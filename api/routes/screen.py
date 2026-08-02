from __future__ import annotations
import io
import numpy as np
import pandas as pd
from fastapi import APIRouter, HTTPException, File, UploadFile, Form
from scipy import stats

router = APIRouter()

@router.post("/screen/analyze")
async def analyze_screen(
    file: UploadFile = File(...),
    guide_col: str = Form(...),
    gene_col: str = Form(...),
    control_col: str = Form(...),
    treatment_col: str = Form(...),
    min_reads: int = Form(30)
):
    content = await file.read()
    try:
        df = pd.read_csv(io.StringIO(content.decode("utf-8")))
    except Exception as e:
        raise HTTPException(400, f"Could not parse CSV: {e}")

    required = {guide_col, gene_col, control_col, treatment_col}
    missing = required - set(df.columns)
    if missing:
        raise HTTPException(400,
            f"Missing columns: {', '.join(sorted(missing))}. "
            f"Available: {', '.join(df.columns)}")

    df = df.rename(columns={
        guide_col: "guide",
        gene_col: "gene",
        control_col: "ctrl",
        treatment_col: "trt"
    })

    df["ctrl"] = pd.to_numeric(df["ctrl"], errors="coerce").fillna(0)
    df["trt"]  = pd.to_numeric(df["trt"],  errors="coerce").fillna(0)
    df = df.dropna(subset=["guide", "gene"])

    n_before  = len(df)
    df = df[df["ctrl"] >= min_reads]
    n_filtered = n_before - len(df)

    if len(df) == 0:
        raise HTTPException(400, "No guides passed the minimum read filter")

    ctrl_total = df["ctrl"].sum() + 1
    trt_total  = df["trt"].sum()  + 1
    df["ctrl_rpm"] = (df["ctrl"] + 1) / ctrl_total * 1e6
    df["trt_rpm"]  = (df["trt"]  + 1) / trt_total  * 1e6
    df["log2fc"]   = np.log2(df["trt_rpm"]) - np.log2(df["ctrl_rpm"])

    mu  = df["log2fc"].mean()
    sig = df["log2fc"].std()
    df["z_score"] = (df["log2fc"] - mu) / sig if sig > 0 else 0.0

    gene_rows = []
    for gene_name, grp in df.groupby("gene"):
        vals = grp["log2fc"].values.astype(float)
        n    = len(vals)
        mean_lfc   = float(np.mean(vals))
        median_lfc = float(np.median(vals))

        if n >= 2:
            _, p = stats.ttest_1samp(vals, 0.0)
            p = float(p)
        else:
            z = (vals[0] - mu) / sig if sig > 0 else 0.0
            p = float(stats.norm.sf(abs(z)) * 2)

        gene_rows.append({
            "gene":          gene_name,
            "n_guides":      int(n),
            "mean_log2fc":   round(mean_lfc,   4),
            "median_log2fc": round(median_lfc, 4),
            "p_value":       round(p, 8),
        })

    gene_df = pd.DataFrame(gene_rows)

    p_arr = gene_df["p_value"].values.astype(float)
    n_g   = len(p_arr)
    order = np.argsort(p_arr)
    adj   = p_arr.copy()
    for rank, idx in enumerate(order):
        adj[idx] = min(1.0, p_arr[idx] * n_g / (rank + 1))
    for i in range(len(order) - 2, -1, -1):
        adj[order[i]] = min(float(adj[order[i]]),
                            float(adj[order[i + 1]]))
    gene_df["fdr"] = np.clip(adj, 0.0, 1.0).round(8)

    def classify(row):
        sig_flag = row["fdr"] < 0.25 or row["p_value"] < 0.05
        if row["mean_log2fc"] < -1.0 and sig_flag:
            return "depleted"
        if row["mean_log2fc"] >  1.0 and sig_flag:
            return "enriched"
        return "neutral"

    gene_df["classification"] = gene_df.apply(classify, axis=1)
    gene_df = gene_df.sort_values("mean_log2fc").reset_index(drop=True)
    gene_df["rank"] = gene_df.index + 1

    depleted_n = int((gene_df["classification"] == "depleted").sum())
    enriched_n = int((gene_df["classification"] == "enriched").sum())
    neutral_n  = int((gene_df["classification"] == "neutral").sum())

    print(f"[SCREEN] genes={n_g} "
          f"p_min={p_arr.min():.4f} p_max={p_arr.max():.4f} "
          f"fdr_min={gene_df['fdr'].min():.4f} "
          f"depleted={depleted_n} enriched={enriched_n}")

    sorted_ctrl = np.sort(df["ctrl"].values.astype(float))
    n_c = len(sorted_ctrl)
    gini = (
        (2 * np.sum((np.arange(1, n_c + 1)) * sorted_ctrl)
         / (n_c * sorted_ctrl.sum())) - (n_c + 1) / n_c
    ) if sorted_ctrl.sum() > 0 else 0.0

    corr = (float(np.corrcoef(df["ctrl_rpm"], df["trt_rpm"])[0, 1])
            if len(df) > 1 else 0.0)

    out_guides = df.rename(columns={
        "ctrl": "control_count",
        "trt":  "treatment_count",
        "ctrl_rpm": "control_rpm",
        "trt_rpm":  "treatment_rpm",
    })[[
        "guide", "gene",
        "control_count", "treatment_count",
        "control_rpm", "treatment_rpm",
        "log2fc", "z_score"
    ]].round(3).to_dict("records")

    return {
        "summary": {
            "total_guides":   int(len(df)),
            "total_genes":    int(n_g),
            "guides_filtered": int(n_filtered),
            "depleted_genes": depleted_n,
            "enriched_genes": enriched_n,
            "neutral_genes":  neutral_n,
            "mean_log2fc":    round(float(df["log2fc"].mean()), 4),
            "std_log2fc":     round(float(df["log2fc"].std()),  4),
        },
        "genes":  gene_df.to_dict("records"),
        "guides": out_guides,
        "qc_metrics": {
            "gini_coefficient":  round(float(gini), 4),
            "correlation":       round(corr, 4),
            "zero_count_guides": int((df["ctrl"] == 0).sum()),
        },
    }

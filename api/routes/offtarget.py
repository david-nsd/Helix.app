from __future__ import annotations

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from helix_core.crisprutils import reverse_complement, sanitize_sequence, validate_sequence

router = APIRouter()


class OffTargetRequest(BaseModel):
    guide: str
    sequence: str
    max_mismatches: int = 4


class OffTargetSite(BaseModel):
    position: int
    sequence: str
    pam: str
    strand: str
    total_mismatches: int
    seed_mismatches: int
    distal_mismatches: int
    risk_score: float
    risk_level: str
    mismatch_map: list[bool]


class OffTargetResponse(BaseModel):
    guide: str
    total_sites: int
    on_target: int
    high_risk: int
    medium_risk: int
    low_risk: int
    sites: list[OffTargetSite]


def _scan_strand(
    guide: str, strand_seq: str, strand: str, max_mm: int
) -> list[OffTargetSite]:
    """Scan a single strand for off-target sites with NGG PAM.

    Args:
        guide: 20nt query guide sequence.
        strand_seq: DNA strand to scan.
        strand: "+" or "-" label for reporting.
        max_mm: Maximum allowed mismatches (total).

    Returns:
        List of OffTargetSite objects in discovery order.
    """
    guide      = guide.upper()
    strand_seq = strand_seq.upper()
    sites: list[OffTargetSite] = []
    n = len(strand_seq)

    for i in range(n - 22):
        window = strand_seq[i:i + 20]
        pam    = strand_seq[i + 20:i + 23]

        if len(pam) < 3 or pam[1] != 'G' or pam[2] != 'G':
            continue

        # Seed = positions 8-19 (12 nt nearest PAM); distal = 0-7 (far from PAM)
        distal_mm = sum(1 for j in range(8)      if window[j] != guide[j])
        seed_mm   = sum(1 for j in range(8, 20)  if window[j] != guide[j])
        total_mm  = seed_mm + distal_mm

        if total_mm > max_mm:
            continue

        mismatch_map = [window[j] != guide[j] for j in range(20)]
        risk_score   = max(0.0, min(1.0, 1.0 - (seed_mm * 0.3 + distal_mm * 0.1)))

        if total_mm == 0:
            risk_level = "on-target"
        elif risk_score >= 0.7:
            risk_level = "high"
        elif risk_score >= 0.4:
            risk_level = "medium"
        else:
            risk_level = "low"

        sites.append(OffTargetSite(
            position=i,
            sequence=window,
            pam=pam,
            strand=strand,
            total_mismatches=total_mm,
            seed_mismatches=seed_mm,
            distal_mismatches=distal_mm,
            risk_score=round(risk_score, 4),
            risk_level=risk_level,
            mismatch_map=mismatch_map,
        ))

    return sites


@router.post("/offtarget")
def find_off_targets(req: OffTargetRequest) -> OffTargetResponse:
    """Find off-target sites for a guide across a genomic sequence.

    Args:
        req: OffTargetRequest with guide, sequence, and max_mismatches.

    Returns:
        OffTargetResponse with sites sorted by risk score (descending).

    Raises:
        HTTPException: 400 if guide is not exactly 20nt or sequence is invalid.
        HTTPException: 500 if scanning fails unexpectedly.
    """
    guide = sanitize_sequence(req.guide.upper())[:20]
    if len(guide) != 20:
        raise HTTPException(
            status_code=400,
            detail=f"Guide must be exactly 20 nt after sanitization, got {len(guide)} nt",
        )

    seq, err = validate_sequence(req.sequence)
    if err:
        raise HTTPException(status_code=400, detail=err)

    print(f"Guide: {guide}")
    print(f"Sequence length: {len(seq)}")

    on_target_found = False
    for i in range(len(seq) - 22):
        window = seq[i:i + 23]
        if not all(c in 'ACGT' for c in window.upper()):
            continue
        pam = window[20:23].upper()
        if pam[1] == 'G' and pam[2] == 'G':
            candidate = window[:20].upper()
            mm = sum(a != b for a, b in zip(guide.upper(), candidate))
            if mm == 0:
                on_target_found = True
                print(f"On-target found at position {i}")
    if not on_target_found:
        print("WARNING: Guide not found in sequence — wrong sequence provided")

    print(f"Max mismatches: {req.max_mismatches}")

    try:
        seq_upper = seq.upper()
        pam_count = sum(
            1 for i in range(len(seq_upper) - 22)
            if seq_upper[i + 21] == 'G' and seq_upper[i + 22] == 'G'
        )
        print(f"PAM sites found: {pam_count}")

        fwd_sites = _scan_strand(guide, seq, "+", req.max_mismatches)
        rev_sites = _scan_strand(guide, reverse_complement(seq), "-", req.max_mismatches)
        all_sites = sorted(fwd_sites + rev_sites, key=lambda s: s.risk_score, reverse=True)
        print(f"Sites passing mismatch filter: {len(all_sites)}")

        return OffTargetResponse(
            guide=guide,
            total_sites=len(all_sites),
            on_target=sum(1 for s in all_sites if s.risk_level == "on-target"),
            high_risk=sum(1 for s in all_sites if s.risk_level == "high"),
            medium_risk=sum(1 for s in all_sites if s.risk_level == "medium"),
            low_risk=sum(1 for s in all_sites if s.risk_level == "low"),
            sites=all_sites,
        )
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Off-target scan failed: {exc}")

# Helix — Project Recap

**Version**: 0.4.0 | **Stack**: Python FastAPI + React 19 + SQLite  
**Purpose**: Full-featured CRISPR bioinformatics suite with ML guide scoring, AI assistant, and experimental design tools.

---

## Project Structure

```
Helix.app/
├── api/                    # FastAPI backend
│   ├── main.py             # App entry point, mounts all routers
│   └── routes/             # 14 endpoint modules
├── frontend/               # React SPA (Vite + Tailwind + Three.js)
│   └── src/
│       ├── pages/          # 14 page components
│       ├── App.jsx         # Router, nav, loading screen
│       ├── api.js          # Axios client wrappers
│       └── store.jsx       # React Context (global state)
├── helix_core/             # 31-module reusable genetics library
├── mobile_api.py           # Separate FastAPI app for mobile clients
├── app.py                  # Legacy Streamlit desktop app (deprecated)
├── helix.db                # SQLite cache + session database
├── requirements.txt
└── Procfile / pyproject.toml
```

---

## Backend API Routes

All endpoints are prefixed `/api/`.

### AI (`api/routes/ai.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/ai/chat` | Claude API chat with CRISPR analysis context injected into system prompt |

**Input**: `{ message, context: { sequence, enzyme, grna_count, top_guide, top_score } }`  
**Output**: `{ reply: string }`  
**Model**: `claude-opus-4-6`, max_tokens=1024

---

### gRNA Analysis (`api/routes/grnas.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/grnas` | Finds and ranks gRNA sites in a sequence |

**Input**: `{ sequence, enzyme, max_guides? }`  
**Output**: `{ grnas: [{ sequence, score, gc, position, strand, risk }], gc_track: [...] }`  
Supports: SpCas9, SaCas9, Cas12a, Cpf1. Scores via ML model with GC fallback.

---

### Scoring (`api/routes/scoring.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/score` | Batch ML scoring of gRNA sequences |

**Input**: `{ guides: [string] }`  
**Output**: `{ scores: [{ guide, score, risk }] }`

---

### Base Editor (`api/routes/base_editor.py`, prefix `/api/baseedit`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/analyze` | Scans for base editing opportunities (CBE/ABE editors) |
| POST | `/prime_editor` | Designs prime editing guides (pegRNA + PE3 nicking guide) |

**Editors supported**: BE3, BE4max, ABE7, ABE8e, NG-CBE  
**Prime editor output**: spacer, PAM, RT template, PBS, full pegRNA, edit preview

---

### Off-Target (`api/routes/offtarget.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/offtarget` | Discovers off-target sites by scanning for PAM matches with mismatches |

**Input**: `{ guide, sequence, max_mismatches? }`  
**Output**: `{ sites: [{ sequence, position, mismatches, mismatch_map, risk }] }`  
Seed region (pos 8–19) weighted higher than distal region (0–7).

---

### Outcome Simulation (`api/routes/outcome.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/outcome/simulate` | Monte Carlo NHEJ/HDR repair outcome simulation |

**Input**: `{ sequence, cut_position, n_simulations, cell_type, has_donor }`  
**Output**: `{ frameshift_pct, nhej_pct, hdr_pct, indel_distribution, microhomology }`

---

### Primers (`api/routes/primers.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/primers/design` | Designs flanking PCR primer pairs and HDR donor template |

**Input**: `{ sequence, cut_position, enzyme }`  
**Output**: `{ forward, reverse, amplicon_size, tm, gc, hdr_template }`

---

### ORFs (`api/routes/orfs.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/orfs` | Finds open reading frames across all 6 reading frames |

**Input**: `{ sequence, min_length? }`  
**Output**: `{ orfs: [{ frame, position, length, protein_preview, protein }] }`

---

### Variants (`api/routes/variants.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/variants` | Compares two sequences and classifies SNP/INS/DEL impact on CRISPR sites |

**Input**: `{ reference, query }`  
**Output**: `{ variants: [{ type, position, impact: disrupts_pam|disrupts_seed|safe }] }`

---

### Protein (`api/routes/protein.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/protein/translate` | Translates DNA to protein in specified frame |
| POST | `/protein/analyze` | Computes MW, pI, hydrophobicity, AA composition |

---

### Genes (`api/routes/genes.py`)
| Method | Path | What it does |
|--------|------|-------------|
| GET | `/genes/common` | Returns hardcoded list of 30+ clinically important human genes |
| GET | `/genes/search` | Queries NCBI for genes by name/organism (7-day cache) |
| GET | `/genes/fetch/{accession}` | Fetches sequence from NCBI GenBank (cached) |

---

### Sequences (`api/routes/sequences.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/sequences/save` | Save a named sequence to helix.db |
| GET | `/sequences` | List all saved sequences |
| GET | `/sequences/{id}` | Fetch one sequence |
| DELETE | `/sequences/{id}` | Delete a sequence |

---

### Sessions (`api/routes/session.py`)
| Method | Path | What it does |
|--------|------|-------------|
| GET | `/sessions` | List all sessions (optional username filter) |
| POST | `/sessions/save` | Save full analysis session JSON to helix.db |
| GET | `/sessions/{id}` | Load a session by ID |

---

### Export (`api/routes/export.py`)
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/export/protocol` | Generates Markdown lab protocol with gRNA table and ordering info |
| POST | `/export/idt_order` | Generates CSV for IDT synthesis order |
| POST | `/export/genbank` | Generates GenBank flat file with annotated gRNA binding sites |

---

### PDB (`api/routes/pdb.py`, async)
| Method | Path | What it does |
|--------|------|-------------|
| GET | `/pdb/search` | Full-text search against RCSB PDB |
| GET | `/pdb/sequence-search` | Sequence similarity search (70% identity, 1e-1 E-value) |
| GET | `/pdb/fetch/{pdb_id}` | Downloads raw PDB file |

---

### Screen Analyzer (`api/routes/screen.py`) ⚠️ BROKEN
| Method | Path | What it does |
|--------|------|-------------|
| POST | `/screen/analyze` | MAGeCK-inspired pooled CRISPR screen statistical analysis |

**Input**: CSV with columns: `guide, gene, control_count, treatment_count`  
**Output**: summary, ranked genes, guide stats, QC metrics (Gini, correlation, FDR)  
**Status**: **BROKEN** — syntax error on line 175 (`cla` incomplete token) prevents import

---

## Frontend Pages

All pages in `frontend/src/pages/`.

| Page | Purpose | Key Features |
|------|---------|-------------|
| **Sandbox.jsx** | Main gRNA analysis interface | Enzyme selector, sequence input, gRNA ranking table, GC track, 3D helix (Three.js), ORF/Variant tabs |
| **BaseEditor.jsx** | Base editor guide design | CBE/ABE editor picker, editing window visualization, efficiency/specificity scores |
| **OutcomeSimulator.jsx** | NHEJ/HDR outcome prediction | Monte Carlo results histogram, frameshift %, microhomology table |
| **OffTarget.jsx** | Off-target site discovery | Guide input, mismatch display, seed/distal risk visualization |
| **PrimerDesigner.jsx** | PCR primer design | Flanking primer pairs, Tm/GC display, HDR donor preview |
| **Orfs.jsx** | ORF finder | 6-frame display, start/stop positions, protein preview |
| **Variants.jsx** | Sequence variant analysis | SNP/INS/DEL comparison, PAM/seed disruption badges |
| **Protein.jsx** | Protein property analysis | DNA → protein translation, MW, pI, hydrophobicity, AA composition chart |
| **Library.jsx** | Gene sequence library | NCBI gene search, common genes list, save/load sequences |
| **ExportCenter.jsx** | Protocol & data export | Download lab protocol PDF, IDT order CSV, GenBank format |
| **AiAssistant.jsx** | Claude AI chat | Context-aware CRISPR Q&A using current analysis state |
| **ScreenAnalyzer.jsx** | Pooled screen analysis | CSV upload, volcano plot (SVG), gene ranking table, QC metrics |
| **Animation.jsx** | 3D helix visualization | Educational animated DNA helix |
| **CrisprGame.jsx** | Educational game | Interactive CRISPR gameplay |

---

## Core Science Modules (`helix_core/`)

| Module | Purpose |
|--------|---------|
| **crisprutils.py** | PAM definitions, site finding, reverse complement, GC track, DNA translation, sequence validation |
| **scoring_model.py** | 103-feature ML guide scoring (XGBoost via `scorer.pkl`); falls back to GC fraction |
| **db.py** | SQLite session and cache management (lazy init, CRUD helpers) |
| **offtarget.py** | Seed-and-extend off-target scanning with k-mer indexing |
| **screen_analysis.py** | Pooled screen statistics: RPM normalization, log2FC, Z-scores, Benjamini-Hochberg FDR, Gini |
| **screen_design.py** | Pooled screen guide pool design (broken imports — see Known Issues) |
| **variants.py** | SNP/indel detection and PAM/seed disruption classification |
| **primer.py** | PCR primer optimization: Tm calculation, hairpin detection, GC constraints |
| **codon.py** | Codon optimization utilities |
| **motifs.py** | Sequence motif analysis |
| **msa_utils.py** | Multiple sequence alignment utilities |
| **protein_tools.py** | Protein property calculations |
| **peptidebuilder.py** | 3D peptide structure generation |
| **structure_viewer.py** | PDB structure visualization helpers |
| **visuals.py** | Matplotlib/Plotly visualization utilities |
| **sonify.py** | Audio synthesis from sequences |
| **animations.py** | Animation utilities |
| **constants.py** | Shared project constants |
| **auth.py** | Authentication helpers |
| **io_utils.py** | File I/O utilities |
| **storage.py** | Local storage helpers |
| **progression.py** | User progression tracking |
| **export_state.py** | Session serialization |
| **editor.py** | Base editing support utilities |
| **ai_stub.py / ai_helper.py** | Claude API integration stubs |

---

## ML Model

- **Current version**: v3.0 (trained on Doench 2016 experimental data)
- **Performance**: Pearson r = 0.79 on held-out test set
- **File**: `helix_core/scorer.pkl` (~2.2 MB, XGBoost pickle)
- **Feature vector**: 103 features per 20nt guide
  - GC%, seed GC, melting temperature, self-complementarity
  - Dinucleotide frequencies (all 16 pairs)
  - One-hot position encoding
- **Scoring function**: `helix_core/scoring_model.py` → `score_guide_ml()`
- **Fallback**: GC fraction when model unavailable

---

## Features Built

| Feature | Status |
|---------|--------|
| gRNA site finding (SpCas9, SaCas9, Cas12a, Cpf1) | Working |
| ML guide scoring (XGBoost, Doench 2016) | Working |
| Off-target analysis | Working |
| Base editor guide design (CBE/ABE) | Working |
| Prime editor design (pegRNA + PE3) | Working |
| NHEJ/HDR outcome simulation (Monte Carlo) | Working |
| PCR primer design | Working |
| ORF finder (6 frames) | Working |
| Variant impact assessment | Working |
| Protein translation + property analysis | Working |
| NCBI gene search + sequence fetch | Working |
| Sequence library (save/load) | Working |
| Analysis session save/load | Working |
| Export center (protocol, IDT CSV, GenBank) | Working |
| PDB structure search + download | Working |
| AI assistant (Claude chat with CRISPR context) | Working |
| 3D helix visualization (Three.js) | Working |
| CRISPR educational game | Working |
| Pooled screen analysis (MAGeCK-inspired) | Broken (syntax error in screen.py:175) |

---

## Known Issues

1. **`api/routes/screen.py` line 175** — Syntax error: incomplete token `cla` at end of line.  
   Screen analysis endpoint (`POST /screen/analyze`) will fail to import and crash the API.

2. **`helix_core/screen_design.py` lines 11–16** — Bare imports (`from crisprutils import ...`, `from offtarget import ...`) will fail when called from the API context. Should be `from helix_core.crisprutils import ...`.

3. **CORS** — All origins allowed (`allow_origins=["*"]`) — acceptable for development, must be restricted before production deployment.

4. **No test suite** — No test files found anywhere in the project.

---

## Tech Stack

### Backend (`requirements.txt`)
| Package | Version | Purpose |
|---------|---------|---------|
| fastapi | 0.118 | API framework |
| uvicorn | 0.37 | ASGI server |
| starlette | 0.48 | ASGI toolkit |
| anthropic | 0.97.0 | Claude API client |
| biopython | 1.85 | Bioinformatics (NCBI, sequence parsing) |
| numpy | 2.3.3 | Numerical computing |
| pandas | 2.3.3 | Data analysis |
| scipy | 1.17.1 | Statistical tests |
| scikit-learn | 1.8.0 | ML utilities |
| matplotlib | 3.7+ | Plotting |
| pydantic | 2.11.9 | Data validation |
| python-dotenv | 1.1.1 | Environment variables |
| PyJWT | 2.10.1 | Auth tokens |
| requests | 2.32.5 | HTTP client |

### Frontend (`package.json`)
| Package | Version | Purpose |
|---------|---------|---------|
| react | 19.2.5 | UI framework |
| react-dom | 19.2.5 | DOM renderer |
| react-router-dom | 7.14.2 | Client-side routing |
| three | 0.184 | 3D helix visualization |
| axios | 1.15.2 | HTTP requests to backend |
| vite | 8.0 | Build tool + dev server |
| tailwindcss | 4.2.4 | Utility CSS |

---

## Deployment

- **Backend**: FastAPI app via `api/main.py`, served with uvicorn (port 8000 default)
- **Frontend**: Vite SPA, built to `frontend/dist/`, served statically from FastAPI root
- **Database**: SQLite at `helix.db` in project root (gene cache, sequence cache, sessions)
- **Environment**: Requires `ANTHROPIC_API_KEY` in `.env`
- **Procfile**: Present (configured for cloud deployment, e.g. Railway/Heroku)
- **Mobile API**: `mobile_api.py` runs as a separate FastAPI process

---

*Generated by reading source code — 2026-08-02*

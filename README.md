# Helix

Status: in development. Built with AI-assisted development (Claude Code).

Helix is a Python bioinformatics project with a Streamlit desktop UI, a small FastAPI mobile API, and shared core genomics utilities. Gradient-boosting regressor (103 sequence features) trained to approximate Azimuth scores from the Doench 2016 dataset (Pearson r = 0.79 on a random held-out split). Not validated against wet-lab results.

## Components

- `app.py`: main Streamlit application
- `mobile_api.py`: FastAPI service for RNA tooling
- `helix_core/`: shared logic for sequence analysis and visualization support
- `helix_desktop/`: desktop UI helpers
- `mobile_app/`: mobile-facing package code

## Configuration

Set these environment variables before starting the Streamlit app:

- `HELIX_USER`
- `HELIX_PASS`
- `HELIX_SALT` (optional)

The project can also load these values from `.env`.

## Run

- Web UI: `helix-web`
- Desktop window on PC: `helix-pc`
- API server: `helix-api`
- Windows shortcut in this repo: `Run Helix Desktop.bat`

### 🎯 ML Scoring Model
- Gradient-boosting regressor (103 sequence features) trained to approximate Azimuth scores from the Doench 2016 dataset (Pearson r = 0.79 on a random held-out split). Not validated against wet-lab results.
- Test set performance: R² = 0.61, Pearson r = 0.79
- 82 sequence features: position-specific one-hot encoding, GC content, seed region GC, melting temperature, self-complementarity, dinucleotide frequencies
- Comparable to published sequence-only CRISPR scoring models
- Reference: Doench JG et al. "Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9." Nature Biotechnology 34, 184–191 (2016)

## Citation

If you use Helix in your research, please cite:

**Helix CRISPR Suite**
David Stoica, https://github.com/david-nsd/Helix.app

**gRNA scoring model training data:**
Doench JG, Fusi N, Sullender M, et al. Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature Biotechnology 34, 184–191 (2016).

**Repair outcome prediction:**
Shen MW, Arbab M, Hsu JY, et al. Predictable and programmable DNA deletions using CRISPR/Cas9. Nature Biotechnology 36, 1060–1068 (2018).

# Epigenetic Age Acceleration as a Modifiable Public Health Target

## A Systematic Review and Meta-Analysis of Environmental, Behavioral, and Social Determinants with Development of MEAB-Index

[![DOI](https://img.shields.io/badge/DOI-10.3390%2Fxxxxx-blue)](https://doi.org/10.3390/xxxxx)
[![PROSPERO](https://img.shields.io/badge/PROSPERO-CRD420261371271-green)](https://www.crd.york.ac.uk/PROSPERO/view/CRD420261371271)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-yellow.svg)](https://www.python.org/)

---

## Overview

This repository contains the complete analytical pipeline, processed datasets, and reproducibility materials for the systematic review and meta-analysis published in:

> **Aliberti, S.M.\*, Marigliano, P.G., & Capunzo, M.** (2026). Epigenetic Age Acceleration as a Modifiable Public Health Target: A Systematic Review and Meta-Analysis of Environmental, Behavioral, and Social Determinants with Development of MEAB-Index. *International Journal of Molecular Sciences*, 27(x), xxxx. https://doi.org/10.3390/xxxxx

The study synthesizes **83 studies** (118 exposure–clock associations) across **four analytical pools** (A–D), quantifying the association between modifiable behavioral, environmental, and social determinants and epigenetic age acceleration (EAA) measured via validated DNA methylation clocks (Horvath, PhenoAge, GrimAge, DunedinPACE).

A novel **Modifiable Epigenetic Aging Burden Index (MEAB-Index)** is introduced to estimate the cumulative preventable burden of biological aging.

---

## Key Findings

| Metric | Estimate | 95% CI | Interpretation |
|---|---|---|---|
| Pool A pooled β | +0.310 years | +0.255 to +0.366 | Significant moderate EAA per unit adverse exposure |
| Pool C pooled OR | 1.749 | 1.348 to 2.269 | Significant increased odds of accelerated EAA |
| MEAB-composite | +0.153 years | +0.077 to +0.230 | Comprehensive burden across all 7 categories |
| MEAB-positive | +0.263 years | +0.158 to +0.368 | Conservative estimate (3 significant categories) |
| CPB | +1.566 years | +1.011 to +2.123 | Cumulative Preventable Burden (upper-bound) |

---

## Repository Structure

```
eaa-meta-analysis/
│
├── README.md                        # This file
├── LICENSE                          # CC BY 4.0
├── CITATION.cff                     # Citation metadata
├── requirements.txt                 # Python dependencies
├── .gitignore
│
├── data/
│   ├── raw/
│   │   ├── extraction_raw_112studies_clean.xlsx
│   │   └── search_records_309.xlsx
│   ├── processed/
│   │   ├── STEP6_standardized.xlsx       # 6 sheets (all_112, pool_A–D, audit)
│   │   ├── pool_A_60studies.csv
│   │   ├── pool_B_9studies.csv
│   │   ├── pool_C_10studies.csv
│   │   ├── pool_D_4studies.csv
│   │   └── excluded_w_rationale.csv
│   └── supplementary/
│       └── Table_S1–S4.xlsx
│
├── scripts/
│   │
│   │  # ── Phase 1: Literature search & screening ──
│   ├── step1a_pubmed_screening.py        # PubMed CSV → screening Excel
│   ├── step1b_scopus_screening.py        # Scopus CSV → screening Excel
│   ├── step2a_autoscreen_pubmed.py       # NLP-assisted triage (PubMed)
│   ├── step2b_autoscreen_scopus.py       # NLP-assisted triage (Scopus)
│   │
│   │  # ── Phase 2: Data extraction & harmonization ──
│   ├── step4_extraction_gpt.py           # PDF→GPT→JSON→Excel extraction
│   ├── step5a_raw_records.py             # Raw record table generation
│   ├── step5b_raw_records_clean.py       # Cleaning, normalization, SE derivation
│   ├── step6_standardization.py          # Effect size harmonization (yi, sei)
│   │
│   │  # ── Phase 3: Quantitative synthesis ──
│   ├── step7_meta_analysis.py            # DerSimonian–Laird model + forest plots
│   ├── step8_bias_sensitivity.py         # Funnel, Egger, Begg, trim-and-fill, LOO
│   └── step9_meab_index.py              # MEAB-Index + bootstrap + sensitivity
│
├── figures/                             # Publication-quality figures (PDF)
├── docs/                                # PRISMA checklist, search strategy, protocol
└── output/                              # Regenerable script outputs
```

---

## Analytical Pipeline

### Phase 1 — Literature Search and Screening

| Step | Script | Description |
|---|---|---|
| 1a | `step1a_pubmed_screening.py` | Imports PubMed CSV (n=159); generates keyword-flagged screening Excel with `decision`, `reason_exclusion` columns |
| 1b | `step1b_scopus_screening.py` | Imports Scopus CSV (n=150); robust record ID (PubMed ID → EID fallback); same screening structure |
| 2a | `step2a_autoscreen_pubmed.py` | NLP-assisted triage: parses abstracts from PubMed export, combines title+abstract keywords, generates `auto_suggestion` (Include/Exclude/Unclear) with confidence. Decision support only — no automatic exclusions |
| 2b | `step2b_autoscreen_scopus.py` | Same NLP triage for Scopus records with abstract from CSV |

### Phase 2 — Data Extraction and Harmonization

| Step | Script | Description |
|---|---|---|
| 4 | `step4_extraction_gpt.py` | PDF → GPT-4.1-mini → structured JSON → append to `extraction_master.xlsx`. Batch-capable (single PDF or folder). Schema: 26 fields per study |
| 5a | `step5a_raw_records.py` | Generates raw extraction table (114 → 112 studies after verification) with summary statistics |
| 5b | `step5b_raw_records_clean.py` | Normalizes study design, country labels; derives SE from CI where missing (`SE = (CI_upper - CI_lower) / 3.92`); generates missing data report |
| 6 | `step6_standardization.py` | **Core harmonization step.** Exhaustive effect-type mapping (50+ unique strings → 5 classes: BETA_UNSTD, BETA_STD, OR, HR, NEEDS_REVIEW). Direction/outcome flags. Computes yi, sei. Assigns to 4 pools. Output: 7-sheet XLSX with full audit trail |

### Phase 3 — Quantitative Synthesis

| Step | Script | Description |
|---|---|---|
| 7 | `step7_meta_analysis.py` | DerSimonian–Laird random-effects model on Pool A (n=60). Overall pooled β, subgroup analyses by 7 exposure categories (k≥3). Forest plots (paginated overall + per-category + summary diamond plot). Excel output with weights |
| 8 | `step8_bias_sensitivity.py` | Publication bias: funnel plots (standard + trim-and-fill), Egger's test, Begg's rank correlation. Sensitivity: leave-one-out analysis (overall + by category). All on Pool A |
| 9 | `step9_meab_index.py` | MEAB-Index computation: composite, positive, CPB. 10,000 bootstrap iterations (percentile CI). Three publication figures: component chart, bootstrap distributions, leave-one-category-out sensitivity. Excel with full MEAB results |

---

## Reproduction Instructions

### Prerequisites

```bash
python --version   # 3.11+

python -m venv venv
source venv/bin/activate

pip install -r requirements.txt
```

### Running the Quantitative Synthesis (Steps 7–9)

These steps require the processed data from Step 6 (`STEP6_standardized.xlsx`):

```bash
cd scripts/
python step7_meta_analysis.py
python step8_bias_sensitivity.py
python step9_meab_index.py
```

Outputs are saved to `step7_outputs/`, `step8_outputs/`, `step9_outputs/` respectively.

### Running the Full Pipeline (Steps 1–9)

Steps 1–5 require the raw search exports and PDF full texts. Step 4 requires an OpenAI API key for GPT-assisted extraction.

```bash
# Phase 1: screening
python step1a_pubmed_screening.py
python step1b_scopus_screening.py
python step2a_autoscreen_pubmed.py
python step2b_autoscreen_scopus.py

# Phase 2: extraction + harmonization (after manual full-text review)
python step4_extraction_gpt.py
python step5a_raw_records.py
python step5b_raw_records_clean.py
python step6_standardization.py

# Phase 3: quantitative synthesis
python step7_meta_analysis.py
python step8_bias_sensitivity.py
python step9_meab_index.py
```

---

## Methodological Notes

### Effect Size Harmonization (Step 6)

Studies were classified into four mutually exclusive analytical pools:

- **Pool A** (primary): Unstandardized β coefficients or mean differences in years of EAA (n=60)
- **Pool B** (sensitivity): Standardized β or SMD in SD units (n=9)
- **Pool C** (secondary): Odds ratios, log-transformed (n=10)
- **Pool D** (secondary): Hazard ratios, log-transformed (n=4)

Standard errors were reconstructed as: `SE = (CI_upper − CI_lower) / 3.92` when not directly reported. For ratio metrics (OR, HR), log transformation was applied before pooling.

### DerSimonian–Laird Estimator

All meta-analytic pooling used the DerSimonian–Laird random-effects estimator, implemented from first principles in Python for full transparency:

```
τ² = max(0, (Q − df) / C)
wᵢ = 1 / (vᵢ + τ²)
θ = Σ(wᵢ · yᵢ) / Σ(wᵢ)
```

### MEAB-Index

Three complementary metrics quantify the cumulative preventable burden:

- **MEAB-composite**: Inverse-variance weighted average across all 7 categories
- **MEAB-positive**: Restricted to statistically significant positive categories
- **CPB**: Unweighted sum of significant positive pooled estimates (upper-bound)

All validated with 10,000 bootstrap iterations (non-parametric, percentile method).

---

## Software and Versions

| Package | Version | Purpose |
|---|---|---|
| Python | ≥ 3.11 | Core language |
| numpy | ≥ 1.24 | Numerical computation |
| pandas | ≥ 2.0 | Data manipulation |
| scipy | ≥ 1.11 | Statistical functions |
| matplotlib | ≥ 3.8 | Publication-quality figures |
| openpyxl | ≥ 3.1 | Excel I/O |
| PyMuPDF (fitz) | ≥ 1.23 | PDF text extraction (Step 4) |
| openai | ≥ 1.0 | GPT-assisted extraction (Step 4) |

---

## Citation

```bibtex
@article{aliberti2026eaa,
  title     = {Epigenetic Age Acceleration as a Modifiable Public Health Target:
               A Systematic Review and Meta-Analysis of Environmental, Behavioral,
               and Social Determinants with Development of {MEAB} Index},
  author    = {Aliberti, Silvana Mirella and Marigliano, Piergiorgio and Capunzo, Mario},
  journal   = {International Journal of Molecular Sciences},
  volume    = {27},
  number    = {x},
  pages     = {xxxx},
  year      = {2026},
  doi       = {10.3390/xxxxx},
  publisher = {MDPI}
}
```

---

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

---

## Contact

- **Silvana Mirella Aliberti** (Corresponding Author) — sialiberti@unisa.it
- **Piergiorgio Marigliano** — piergiorgio.marigliano@uniecampus.it
- **Mario Capunzo** — mcapunzo@unisa.it

Department of Medicine, Surgery and Dentistry "Scuola Medica Salernitana", University of Salerno, Italy

# Data Directory

This directory contains all datasets produced during the systematic review and meta-analysis, organized by processing stage.

## `raw/` — Pre-standardization data

### `01_search_exports/` — Original database exports (April 7, 2026)
| File | Description | Records |
|------|-------------|---------|
| `pubmed_export_159records.csv` | PubMed/MEDLINE raw CSV export | 159 |
| `scopus_export_150records.csv` | Scopus raw CSV export | 150 |

### `02_screening/` — Title/abstract screening and deduplication
| File | Description | Records |
|------|-------------|---------|
| `screening_pubmed_manual.xlsx` | PubMed screening with decision columns | 159 |
| `screening_pubmed_159_master.xlsx` | PubMed master with titles | 159 |
| `screening_scopus_manual.xlsx` | Scopus screening with decision columns | 150 |
| `screening_scopus_autosuggest.xlsx` | Scopus NLP-assisted triage | 150 |
| `screening_combined_303_dedup.xlsx` | Combined PubMed+Scopus after deduplication | 303 |
| `fulltext_eligible_203records.xlsx` | Records advancing to full-text assessment | 203 |

### `03_extraction/` — Data extraction pipeline
| File | Description | Records |
|------|-------------|---------|
| `extraction_master.xlsx` | GPT-assisted extraction master (raw) | ~114 |
| `extraction_master.csv` | CSV version of extraction master | ~114 |
| `extraction_master_clean.csv` | Cleaned extraction master | ~114 |
| `batch_review_40studies.xlsx` | Batch verification subset (40 studies) | 40 |
| `extraction_raw_114studies.xlsx` | Raw extraction (114 studies) | 114 |
| `extraction_raw_114studies_clean.xlsx` | Cleaned with normalized design/country | 114 |
| `extraction_raw_112studies_clean.xlsx` | Final after residual duplicate removal | 112 |
| `extraction_missing_report.xlsx` | Missing data report (SE, CI, p-value) | — |

## `processed/` — Analysis-ready data (post-Step 6)

| File | Description | Details |
|------|-------------|---------|
| `STEP6_standardized.xlsx` | **Master standardized file** | 6 sheets: `all_112_final` (112 studies, 39 vars), `pool_A_primary_beta` (60), `pool_B_std_beta` (9), `pool_C_OR_log` (10), `pool_D_HR_log` (4), `excluded_w_rationale` (29) |
| `pool_A_60studies.csv` | Pool A extract (CSV) | 60 studies, unstandardized β, years EAA |
| `pool_B_9studies.csv` | Pool B extract (CSV) | 9 studies, standardized β/SMD |
| `pool_C_10studies.csv` | Pool C extract (CSV) | 10 studies, log(OR) |
| `pool_D_4studies.csv` | Pool D extract (CSV) | 4 studies, log(HR) |
| `excluded_29studies_w_rationale.csv` | Excluded from quantitative synthesis | 29 studies with documented rationale |
| `step7_meta_analysis_results.xlsx` | Step 7 output | Pooled estimates, study weights, narrative-only studies |
| `step8_bias_sensitivity_results.xlsx` | Step 8 output | Egger/Begg tests, trim-and-fill, leave-one-out |

## Data Provenance

All data are derived from publicly available published sources cited in the manuscript. The complete dataset is available from the corresponding author upon reasonable request.

## Key Variables in Pool Files

| Variable | Description |
|----------|-------------|
| `study_id` | Unique study identifier |
| `first_author` | First author surname |
| `year` | Publication year |
| `exposure_category` | Harmonized exposure domain |
| `clock_type` | Epigenetic clock used |
| `yi` | Standardized effect size (common metric per pool) |
| `sei` | Standard error of yi |
| `effect_class` | Pool assignment (BETA_UNSTD, BETA_STD, OR, HR) |
| `flag_eaa_as_exposure` | True if EAA used as predictor (wrong direction) |
| `flag_wrong_outcome` | True if outcome is not EAA |
| `yi_sei_formula` | Audit trail: formula used for yi/sei computation |

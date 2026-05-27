# Analytical Protocol

## Pre-specified Analytical Framework

This document describes the analytical protocol defined prior to quantitative data synthesis, as registered in PROSPERO (CRD420261371271).

### Primary Research Question

What is the magnitude and direction of the association between modifiable behavioral, environmental, and social determinants and epigenetic age acceleration (EAA) measured using validated DNA methylation-based clocks in adult populations?

### Eligibility Criteria (PICOS)

- **Population**: Adults (≥18 years) from general population or community-based samples
- **Exposure**: Modifiable behavioral, environmental, psychosocial, or socioeconomic determinants
- **Comparator**: Different levels of exposure or reference/unexposed groups
- **Outcome**: EAA measured using validated clocks (Horvath, Hannum, PhenoAge, GrimAge/GrimAge2, DunedinPoAm, DunedinPACE)
- **Study design**: Observational (cross-sectional, cohort, case-control) and interventional

### Effect Size Harmonization Framework

Four mutually exclusive analytical pools based on statistical format:

| Pool | Metric | Scale | Role |
|------|--------|-------|------|
| A | Unstandardized β / mean difference | Years of EAA | Primary analysis |
| B | Standardized β / SMD | SD units | Sensitivity |
| C | Odds Ratio (log-transformed) | log(OR) | Secondary |
| D | Hazard Ratio (log-transformed) | log(HR) | Secondary |

### SE Reconstruction

When standard errors were not reported:
- Linear metrics: `SE = (CI_upper − CI_lower) / 3.92`
- Ratio metrics: `SE_log = (log(CI_upper) − log(CI_lower)) / 3.92`

### Meta-analytic Model

- **Estimator**: DerSimonian–Laird (random-effects)
- **Heterogeneity**: Cochran's Q (p < 0.10), I², τ²
- **Subgroup analyses**: By exposure category (k ≥ 3 required)
- **Meta-regression**: Exposure category, clock type, study design, adjustment level
- **Publication bias**: Funnel plots, Egger's test, Begg's test, trim-and-fill
- **Sensitivity**: Leave-one-out, restriction to fully adjusted models

### MEAB-Index (Novel)

Three complementary metrics derived from Pool A subgroup estimates:

1. **MEAB-composite**: Inverse-variance weighted average across all categories
2. **MEAB-positive**: Restricted to significant positive categories
3. **CPB**: Unweighted sum of significant positive estimates

Uncertainty: 10,000 bootstrap iterations (percentile method).

### Deviations from Registered Protocol

1. Title refined to better reflect final scope
2. MEAB-Index added (not in original registration)

All deviations defined prior to quantitative synthesis.

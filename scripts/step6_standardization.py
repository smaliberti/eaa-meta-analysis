"""
step6_standardization.py
========================
Meta-analysis: Modifiable determinants of epigenetic age acceleration
Step 6: Standardization of effect sizes into common metric (yi, sei)

METHODOLOGICAL DECISIONS (document for Methods section):
---------------------------------------------------------
1. PRIMARY METRIC: Unstandardized beta coefficient (years of EAA per unit exposure)
   - Directly interpretable: "X more/fewer years of biological aging per unit increase in exposure"
   - Applied to: beta coefficients, mean differences, regression coefficients in year scale

2. SECONDARY METRICS (separate pools, not combined with primary):
   - Pool B: Standardized beta (SD units) – different scale, sensitivity analysis
   - Pool C: log(OR) – EAA dichotomized as outcome, logistic regression studies
   - Pool D: log(HR) – time-to-event analyses

3. DIRECTION FLAG: Studies where EAA is the PREDICTOR (not outcome) are flagged and
   excluded from primary analysis. These have "epigenetic age acceleration" as exposure_category
   and a downstream health outcome (mortality, dementia, etc.) – reversed direction.

4. OUTCOME FLAG: Studies where outcome is NOT EAA (e.g., dementia cases, cognitive function,
   mortality) are flagged and excluded from primary analysis.

5. SE RECONSTRUCTION: SE = (CI_upper - CI_lower) / 3.92 for ratio metrics on original scale;
   SE_log = (log(CI_upper) - log(CI_lower)) / 3.92 for ratio metrics after log transformation.

6. SPECIAL CASES requiring manual review are flagged with NEEDS_REVIEW status.

Authors: [your name], supervised by [professor name]
Date: 2026
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# SECTION 1: LOAD DATA
# ============================================================

INPUT_FILE  = "extraction_raw_112studies_clean.xlsx"
OUTPUT_FILE = "step6_standardized.xlsx"

df = pd.read_excel(INPUT_FILE)
print(f"Loaded {len(df)} studies from {INPUT_FILE}")
print(f"Columns: {list(df.columns)}\n")

# Work on a copy, preserve all original columns
out = df.copy()


# ============================================================
# SECTION 2: EXHAUSTIVE EFFECT TYPE MAPPING
# Every single unique effect_type is mapped explicitly.
# This ensures no study is silently misclassified.
# ============================================================

# Map from raw effect_type string  →  (effect_class, notes)
# effect_class values:
#   BETA_UNSTD  = unstandardized beta / mean difference (primary pool)
#   BETA_STD    = standardized beta / SMD (sensitivity pool B)
#   OR          = odds ratio (secondary pool C)
#   HR          = hazard ratio (secondary pool D)
#   NEEDS_REVIEW = requires manual adjudication

EFFECT_TYPE_MAP = {

    # ---- BETA_UNSTD (primary) -------------------------------------------
    "Adjusted beta coefficient":
        ("BETA_UNSTD", "Adjusted beta in years; use as-is"),

    "Beta coefficient":
        ("BETA_UNSTD", "Linear regression beta"),

    "Beta coefficient (IVW MR estimate)":
        ("BETA_UNSTD", "Mendelian Randomisation IVW; causal beta estimate; flag MR for sensitivity"),

    "Beta coefficient (difference in years or SD units of biological aging)":
        ("BETA_UNSTD", "Reported in years; SD-unit cases need sub-flag — check effect_value scale"),

    "Beta coefficient (linear regression)":
        ("BETA_UNSTD", "Standard linear beta"),

    "Beta coefficient (years change in age acceleration per interquintile increase in diet score)":
        ("BETA_UNSTD", "Years per quintile increment; consistent with primary metric"),

    "Beta coefficient (years change in epigenetic age per doubling of intake or per unit increase in ratio)":
        ("BETA_UNSTD", "Years per doubling; note exposure unit in conversion_notes"),

    "Beta coefficient (years increase in GrimAge EAA per SD increase in exposure)":
        ("BETA_UNSTD", "Years per SD exposure; exposure unit is SD not raw — note in conversion_notes"),

    "Beta coefficient (years increase in epigenetic age per unit increase in TCDD)":
        ("BETA_UNSTD", "Years per TCDD unit; primary metric"),

    "Change in epigenetic age acceleration (years)":
        ("BETA_UNSTD", "Mean difference / change score in years"),

    "Difference in epigenetic age acceleration (years) or pace of aging":
        ("BETA_UNSTD", "Years difference; DunedinPACE pace-of-aging part is dimensionless — sub-flag if pace"),

    "Difference in epigenetic age acceleration (years) or rate per chronological year":
        ("BETA_UNSTD", "Years difference"),

    "Difference in years of epigenetic age acceleration per 10-unit increase in BDI-II score":
        ("BETA_UNSTD", "Years per 10 BDI-II points; note exposure unit"),

    "Mean difference in DNAm age acceleration between active and inactive co-twins":
        ("BETA_UNSTD", "Twin-design mean difference in years"),

    "Mean difference in age acceleration (years) by social mobility status":
        ("BETA_UNSTD", "Mean difference in years"),

    "Mean difference in epigenetic age (years)":
        ("BETA_UNSTD", "Mean difference in years"),

    "Regression coefficient (B)":
        ("BETA_UNSTD", "Unstandardized regression coefficient"),

    "Regression coefficient (beta)":
        ("BETA_UNSTD", "Linear regression beta"),

    "Regression coefficient (β)":
        ("BETA_UNSTD", "Linear regression beta"),

    "Regression coefficient (β) for linear and quadratic terms":
        ("BETA_UNSTD", "Linear term extracted; quadratic term — verify which was extracted"),

    "Regression coefficient (β) per 50-unit increase in SII":
        ("BETA_UNSTD", "Years per 50-unit SII; note non-standard exposure unit"),

    "Regression coefficient (β) representing years of increase/decrease in epigenetic age acceleration or rate of aging":
        ("BETA_UNSTD", "Years EAA; primary metric"),

    "absolute change in epigenetic age acceleration per IQR increase in air pollutant concentration":
        ("BETA_UNSTD", "Years per IQR; note IQR-scaled exposure"),

    "beta coefficient":
        ("BETA_UNSTD", "Linear beta"),

    "beta coefficient (linear regression)":
        ("BETA_UNSTD", "Linear beta"),

    "beta coefficient (weeks)":
        ("BETA_UNSTD", "CAUTION: unit is weeks not years; conversion_notes set; convert /52 if clock in weeks"),

    "beta coefficient (years for KDM-BA and PhenoAge acceleration; percent change for telomere length)":
        ("BETA_UNSTD", "Years part extracted for PhenoAge/KDM; telomere % part flagged separately"),

    "beta coefficient (years increase in BioAgeAccel per unit increase in lnCMI)":
        ("BETA_UNSTD", "Years per ln-unit CMI; primary metric"),

    "beta coefficient (years increase in EAA per additional drink per day)":
        ("BETA_UNSTD", "Years per drink/day; primary metric"),

    "beta coefficient (years increase in GrimAgeAccel per unit increase in alcohol consumption)":
        ("BETA_UNSTD", "Years per alcohol unit; primary metric"),

    "beta coefficient (years increase in biological age acceleration per 1-point SDH increment)":
        ("BETA_UNSTD", "Years per SDH point; primary metric"),

    "beta coefficient (years of age acceleration)":
        ("BETA_UNSTD", "Years; primary metric"),

    "beta coefficient (years of epigenetic age acceleration)":
        ("BETA_UNSTD", "Years; primary metric"),

    "beta coefficient (years or % change per 1 SD increment in diet score)":
        ("BETA_UNSTD", "Years part extracted; verify in conversion_notes"),

    "beta coefficient (years or unit change in epigenetic age acceleration)":
        ("BETA_UNSTD", "Years EAA; primary metric"),

    "beta coefficient per IQR increase in pollutant":
        ("BETA_UNSTD", "Years per IQR; note IQR-scaled exposure"),

    "change in EAA (years) per interquartile range increase in greenness":
        ("BETA_UNSTD", "Years per IQR; primary metric"),

    "linear regression coefficient (beta)":
        ("BETA_UNSTD", "Linear beta"),

    "mean difference":
        ("BETA_UNSTD", "Mean difference — assume years given EAA context"),

    "mean difference in epigenetic age acceleration (years) or estimate difference (DunedinPACE)":
        ("BETA_UNSTD", "Years / DunedinPACE rate; primary metric"),

    "mean difference in gestational length (days) per SD increase in maternal epigenetic age acceleration":
        ("NEEDS_REVIEW", "WRONG DIRECTION: EAA is exposure; outcome is gestational length — exclude from primary"),

    "mean difference in years (or kb for DNAmTL) of epigenetic age acceleration relative to non-smokers":
        ("BETA_UNSTD", "Years EAA relative to non-smokers; DNAmTL kb part — note in conversion"),

    "mean difference in years of EAA":
        ("BETA_UNSTD", "Mean difference in years; primary metric"),

    "regression coefficient":
        ("BETA_UNSTD", "Unstandardized regression coefficient"),

    "regression coefficient (B)":
        ("BETA_UNSTD", "Unstandardized B"),

    "regression coefficient (B) for mean cognitive performance and intraindividual variability":
        ("NEEDS_REVIEW", "WRONG DIRECTION: outcome appears to be cognitive performance, not EAA — verify"),

    "regression coefficient (b)":
        ("BETA_UNSTD", "Unstandardized b"),

    "regression coefficient (beta)":
        ("BETA_UNSTD", "Linear beta"),

    "regression coefficient (β)":
        ("BETA_UNSTD", "Linear beta"),

    # ---- BETA_STD (sensitivity pool B) ---------------------------------
    "Beta coefficient (standardized effect size)":
        ("BETA_STD", "Standardized beta; SD units; sensitivity pool B"),

    "Standardized beta coefficient":
        ("BETA_STD", "Standardized beta; sensitivity pool B"),

    "Standardized beta coefficient (SD change in epigenetic age acceleration per SD increase in anxiety)":
        ("BETA_STD", "SD-SD standardized; sensitivity pool B"),

    "Standardized beta coefficient (per 1 SD increase in physical activity score)":
        ("BETA_STD", "SD-standardized; sensitivity pool B"),

    "Standardized mean difference":
        ("BETA_STD", "SMD (Cohen's d equivalent); sensitivity pool B"),

    "beta coefficient (standardized)":
        ("BETA_STD", "Standardized beta; sensitivity pool B"),

    "standardized beta":
        ("BETA_STD", "Standardized beta; sensitivity pool B"),

    "standardized beta coefficient (SD units of PhenoAge advancement)":
        ("BETA_STD", "SD units; sensitivity pool B"),

    "standardized regression coefficient":
        ("BETA_STD", "Standardized; sensitivity pool B"),

    # ---- OR (secondary pool C) ------------------------------------------
    "Odds Ratio (OR) for MetS-ATP per unit increase in UA and UAA":
        ("NEEDS_REVIEW", "WRONG DIRECTION: EAA (UA/UAA) is EXPOSURE predicting MetS — exclude from primary"),

    "Odds Ratio (OR) for PhenoAgeAccel":
        ("OR", "OR for dichotomized PhenoAgeAccel; secondary pool C; yi=log(OR)"),

    "Odds ratio":
        ("OR", "Logistic regression OR; secondary pool C; yi=log(OR)"),

    "Odds ratio per 1 SD increase":
        ("OR", "OR per 1 SD exposure; secondary pool C; yi=log(OR)"),

    "Odds ratio per 1 standard deviation increase":
        ("OR", "OR per SD; secondary pool C; yi=log(OR)"),

    "adjusted odds ratio":
        ("OR", "Adjusted OR; secondary pool C; yi=log(OR)"),

    "odds ratio":
        ("OR", "OR; secondary pool C; yi=log(OR)"),

    "odds ratio per 1 SD increase":
        ("OR", "OR per SD; secondary pool C; yi=log(OR)"),

    # ---- HR (secondary pool D) ------------------------------------------
    "Hazard Ratio":
        ("HR", "Cox model HR; secondary pool D; yi=log(HR)"),

    "Hazard Ratio (HR) for mediation indirect effect":
        ("NEEDS_REVIEW", "Mediation indirect HR — not directly poolable; flag for narrative synthesis"),

    "Hazard ratio":
        ("HR", "Cox model HR; secondary pool D; yi=log(HR)"),

    "Hazard ratio per 5-year increase in EAA (10% increase for DunedinPoAm)":
        ("NEEDS_REVIEW", "WRONG DIRECTION: EAA is EXPOSURE predicting outcome (mortality/disease) — exclude"),

    "hazard ratio":
        ("HR", "Cox HR; secondary pool D; yi=log(HR)"),

    "hazard ratio (HR) for mortality; odds ratio (OR) for cancer prevalence":
        ("NEEDS_REVIEW", "WRONG DIRECTION: EAA predicts mortality/cancer — exclude from primary; narrative only"),

    "hazard ratio per 5-year increase in PhenoAgeAccel":
        ("NEEDS_REVIEW", "WRONG DIRECTION: EAA (PhenoAgeAccel) is EXPOSURE — exclude"),

    "hazard ratio per standard deviation increase in epigenetic age acceleration":
        ("NEEDS_REVIEW", "WRONG DIRECTION: EAA is EXPOSURE — exclude"),

    # ---- SPECIAL / NEEDS REVIEW -----------------------------------------
    "Machine learning model variable importance and interaction effects; linear regression interaction term for diabetes and blood toluene":
        ("NEEDS_REVIEW", "ML variable importance is not a standard effect size; linear interaction term may be extractable — check paper"),

    "Regression coefficient or Odds Ratio (OR) for 4-year mortality":
        ("NEEDS_REVIEW", "WRONG DIRECTION: EAA predicts 4-year mortality — exclude"),

    "Standardized beta coefficient; Hazard ratio":
        ("NEEDS_REVIEW", "Mixed report: extract standardized beta (pool B) OR HR (pool D); check which refers to EAA as outcome"),

    "beta coefficient and odds ratio":
        ("BETA_UNSTD", "Extract beta for primary pool; OR part → pool C separately; note in conversion"),

    "beta coefficient (linear regression) and odds ratio (logistic regression) for accelerated aging":
        ("BETA_UNSTD", "Extract linear beta for primary pool; OR part → pool C; note in conversion"),

    "effect size (rbc) for age acceleration difference between regions":
        ("NEEDS_REVIEW", "Rank-biserial correlation (rbc); non-standard — convert to d = 2r/sqrt(1-r²) if needed; review"),

    "mediation effect":
        ("NEEDS_REVIEW", "Indirect (mediation) effect — not directly poolable with direct associations; narrative synthesis"),

    "absolute difference in dementia cases per 10,000 person-years":
        ("NEEDS_REVIEW", "WRONG OUTCOME: outcome is dementia incidence, not EAA — exclude from quantitative synthesis"),

    "standardized beta":
        ("BETA_STD", "Standardized beta; sensitivity pool B"),
}


# ============================================================
# SECTION 3: APPLY MAPPING
# ============================================================

def get_class_and_notes(et):
    """Look up effect_type in map. Return (class, notes) or NEEDS_REVIEW if not found."""
    if pd.isna(et):
        return ("NEEDS_REVIEW", "MISSING effect_type — manual check required")
    et_str = str(et).strip()
    if et_str in EFFECT_TYPE_MAP:
        return EFFECT_TYPE_MAP[et_str]
    # Fuzzy fallback: try lowercase match
    for key, val in EFFECT_TYPE_MAP.items():
        if key.lower() == et_str.lower():
            return val
    return ("NEEDS_REVIEW", f"Unrecognized effect_type: [{et_str}] — add to map after manual review")


classes_notes = out['effect_type'].apply(get_class_and_notes)
out['effect_class']       = classes_notes.apply(lambda x: x[0])
out['conversion_notes']   = classes_notes.apply(lambda x: x[1])


# ============================================================
# SECTION 4: DIRECTION AND OUTCOME FLAGS
# Studies where EAA is the predictor (not the outcome) must
# be excluded from the primary meta-analytic pool.
# These are identifiable by exposure_category containing
# epigenetic age / biological aging terminology.
# ============================================================

EAA_AS_EXPOSURE_PATTERNS = [
    'epigenetic age acceleration',
    'biological aging',
    'biological ageing',
    'phenotypic age acceleration',
    'epigenetic aging',
    'biological age',
    'urological age',
]

def is_eaa_as_exposure(row):
    """Flag if EAA/biological aging is listed as the EXPOSURE (wrong direction)."""
    exp_cat = str(row.get('exposure_category', '')).lower()
    for pat in EAA_AS_EXPOSURE_PATTERNS:
        if pat in exp_cat:
            return True
    # Additional explicit check from NEEDS_REVIEW notes
    notes = str(row.get('conversion_notes', '')).lower()
    if 'wrong direction' in notes:
        return True
    return False

def is_wrong_outcome(row):
    """Flag if the outcome in the paper is NOT EAA (e.g., dementia, mortality, cognition)."""
    notes = str(row.get('conversion_notes', '')).lower()
    return 'wrong outcome' in notes or 'wrong direction' in notes

out['flag_eaa_as_exposure'] = out.apply(is_eaa_as_exposure, axis=1)
out['flag_wrong_outcome']   = out.apply(is_wrong_outcome, axis=1)


# ============================================================
# SECTION 5: COMPUTE yi AND sei
# ============================================================

def compute_yi_sei(row):
    """
    Apply the correct yi/sei formula based on effect_class.

    Returns:
        yi              : standardized effect size (common metric)
        sei             : standard error of yi
        log_transformed : True if log was applied (HR, OR)
        qi_formula      : human-readable formula used (for audit trail)
    """
    ec      = row['effect_class']
    ev      = row['effect_value']
    se_d    = row['SE_derived']       # already computed as (CI_upper-CI_lower)/3.92
    ci_lo   = row['CI_lower']
    ci_hi   = row['CI_upper']
    se_orig = row['SE']               # SE as reported in original paper (rare)

    # Use reported SE if available and credible; else use SE_derived
    def best_se_linear():
        if pd.notna(se_orig) and float(se_orig) > 0:
            return float(se_orig), "SE from paper"
        if pd.notna(se_d) and float(se_d) > 0:
            return float(se_d), "SE derived from CI"
        return np.nan, "SE unavailable"

    yi  = np.nan
    sei = np.nan
    log_t = False
    formula = "not computed"

    try:
        ev_f = float(ev)
    except (TypeError, ValueError):
        # Handles dict-strings (e.g., nihms-1835407 Spartano)
        return np.nan, np.nan, False, "effect_value is non-numeric (dict/string) — NEEDS_REVIEW"

    if ec == 'BETA_UNSTD':
        yi = ev_f
        sei_val, src = best_se_linear()
        sei = sei_val
        log_t = False
        formula = f"yi = beta; sei = {src}"

    elif ec == 'BETA_STD':
        yi = ev_f
        sei_val, src = best_se_linear()
        sei = sei_val
        log_t = False
        formula = f"yi = std_beta (SD units); sei = {src}"

    elif ec == 'OR':
        if ev_f <= 0:
            return np.nan, np.nan, False, "OR <= 0: invalid value"
        yi = np.log(ev_f)
        log_t = True
        if pd.notna(ci_lo) and pd.notna(ci_hi) and float(ci_lo) > 0 and float(ci_hi) > 0:
            sei = (np.log(float(ci_hi)) - np.log(float(ci_lo))) / 3.92
            formula = "yi = log(OR); sei = (log(CI_hi)-log(CI_lo))/3.92"
        elif pd.notna(se_orig) and float(se_orig) > 0:
            sei = float(se_orig)
            formula = "yi = log(OR); sei = SE from paper (log scale assumed)"
        else:
            sei = np.nan
            formula = "yi = log(OR); sei = unavailable"

    elif ec == 'HR':
        if ev_f <= 0:
            return np.nan, np.nan, False, "HR <= 0: invalid value"
        yi = np.log(ev_f)
        log_t = True
        if pd.notna(ci_lo) and pd.notna(ci_hi) and float(ci_lo) > 0 and float(ci_hi) > 0:
            sei = (np.log(float(ci_hi)) - np.log(float(ci_lo))) / 3.92
            formula = "yi = log(HR); sei = (log(CI_hi)-log(CI_lo))/3.92"
        elif pd.notna(se_orig) and float(se_orig) > 0:
            sei = float(se_orig)
            formula = "yi = log(HR); sei = SE from paper (log scale assumed)"
        else:
            sei = np.nan
            formula = "yi = log(HR); sei = unavailable"

    elif ec == 'NEEDS_REVIEW':
        formula = "NOT COMPUTED — requires manual adjudication"

    return yi, sei, log_t, formula


results = out.apply(compute_yi_sei, axis=1, result_type='expand')
results.columns = ['yi', 'sei', 'log_transformed', 'yi_sei_formula']
out = pd.concat([out, results], axis=1)


# ============================================================
# SECTION 6: POOLABLE FLAG
# A study is poolable in the PRIMARY pool if:
#   1. effect_class is BETA_UNSTD (or MD — already mapped there)
#   2. NOT flagged as wrong direction
#   3. NOT flagged as wrong outcome
#   4. yi and sei are both numeric and finite
# ============================================================

def is_poolable_primary(row):
    if row['effect_class'] != 'BETA_UNSTD':
        return False
    if row['flag_eaa_as_exposure']:
        return False
    if row['flag_wrong_outcome']:
        return False
    if pd.isna(row['yi']) or pd.isna(row['sei']):
        return False
    if not np.isfinite(row['yi']) or not np.isfinite(row['sei']):
        return False
    return True

def is_poolable_std(row):
    if row['effect_class'] != 'BETA_STD':
        return False
    if row['flag_eaa_as_exposure'] or row['flag_wrong_outcome']:
        return False
    return pd.notna(row['yi']) and pd.notna(row['sei'])

def is_poolable_or(row):
    if row['effect_class'] != 'OR':
        return False
    if row['flag_eaa_as_exposure'] or row['flag_wrong_outcome']:
        return False
    return pd.notna(row['yi']) and pd.notna(row['sei'])

def is_poolable_hr(row):
    if row['effect_class'] != 'HR':
        return False
    if row['flag_eaa_as_exposure'] or row['flag_wrong_outcome']:
        return False
    return pd.notna(row['yi']) and pd.notna(row['sei'])

out['poolable_primary']   = out.apply(is_poolable_primary, axis=1)
out['poolable_pool_B_std']= out.apply(is_poolable_std, axis=1)
out['poolable_pool_C_OR'] = out.apply(is_poolable_or, axis=1)
out['poolable_pool_D_HR'] = out.apply(is_poolable_hr, axis=1)


# ============================================================
# SECTION 7: PRINT SUMMARY REPORT
# ============================================================

print("=" * 65)
print("STEP 6 — STANDARDIZATION REPORT")
print("=" * 65)

print("\n--- Effect class distribution ---")
print(out['effect_class'].value_counts().to_string())

print("\n--- Direction / outcome flags ---")
print(f"  EAA as exposure (wrong direction): {out['flag_eaa_as_exposure'].sum()}")
print(f"  Wrong outcome (not EAA):           {out['flag_wrong_outcome'].sum()}")

print("\n--- Poolability summary ---")
print(f"  Pool A – PRIMARY (BETA_UNSTD, correct dir.): {out['poolable_primary'].sum()}")
print(f"  Pool B – sensitivity (BETA_STD):             {out['poolable_pool_B_std'].sum()}")
print(f"  Pool C – secondary (OR, log scale):          {out['poolable_pool_C_OR'].sum()}")
print(f"  Pool D – secondary (HR, log scale):          {out['poolable_pool_D_HR'].sum()}")
print(f"  NEEDS_REVIEW:                                {(out['effect_class'] == 'NEEDS_REVIEW').sum()}")

print("\n--- yi / sei completeness in primary pool ---")
prim = out[out['poolable_primary']]
print(f"  Studies in primary pool: {len(prim)}")
print(f"  yi available:  {prim['yi'].notna().sum()}")
print(f"  sei available: {prim['sei'].notna().sum()}")

print("\n--- NEEDS_REVIEW studies (require manual adjudication) ---")
needs = out[out['effect_class'] == 'NEEDS_REVIEW']
for _, r in needs.iterrows():
    print(f"  {r['study_id'][:50]:50s} | {r['first_author']:20s} | {str(r['conversion_notes'])[:80]}")

print("\n--- Studies with MISSING sei in primary pool ---")
missing_sei = out[out['poolable_primary'] & out['sei'].isna()]
if len(missing_sei) == 0:
    print("  None — all primary pool studies have sei ✓")
else:
    for _, r in missing_sei.iterrows():
        print(f"  {r['study_id'][:50]:50s} | sei formula: {r['yi_sei_formula']}")


# ============================================================
# SECTION 8: SAVE OUTPUT
# ============================================================

# Column order: original columns first, then new computed columns
new_cols = [
    'effect_class',
    'flag_eaa_as_exposure',
    'flag_wrong_outcome',
    'yi',
    'sei',
    'log_transformed',
    'poolable_primary',
    'poolable_pool_B_std',
    'poolable_pool_C_OR',
    'poolable_pool_D_HR',
    'yi_sei_formula',
    'conversion_notes',
]

# Put new cols right after 'SE_derived' for readability
orig_cols = list(df.columns)
insert_after = orig_cols.index('SE_derived') + 1
final_cols = orig_cols[:insert_after] + new_cols + orig_cols[insert_after:]
# Drop duplicates if any col appeared twice
seen = set()
final_cols_dedup = []
for c in final_cols:
    if c not in seen:
        final_cols_dedup.append(c)
        seen.add(c)

out_final = out[final_cols_dedup]

with pd.ExcelWriter(OUTPUT_FILE, engine='openpyxl') as writer:

    # Sheet 1: Full standardized dataset
    out_final.to_excel(writer, sheet_name='all_112_standardized', index=False)

    # Sheet 2: Primary pool only (ready for random-effects model)
    primary = out_final[out_final['poolable_primary']].copy()
    primary.to_excel(writer, sheet_name='pool_A_primary_beta', index=False)

    # Sheet 3: Standardized beta pool (sensitivity)
    pool_b = out_final[out_final['poolable_pool_B_std']].copy()
    pool_b.to_excel(writer, sheet_name='pool_B_std_beta', index=False)

    # Sheet 4: OR pool (secondary)
    pool_c = out_final[out_final['poolable_pool_C_OR']].copy()
    pool_c.to_excel(writer, sheet_name='pool_C_OR_log', index=False)

    # Sheet 5: HR pool (secondary)
    pool_d = out_final[out_final['poolable_pool_D_HR']].copy()
    pool_d.to_excel(writer, sheet_name='pool_D_HR_log', index=False)

    # Sheet 6: NEEDS_REVIEW — manual adjudication list
    review = out_final[out_final['effect_class'] == 'NEEDS_REVIEW'].copy()
    review.to_excel(writer, sheet_name='NEEDS_REVIEW_manual', index=False)

    # Sheet 7: Audit log — every study with its class, flags, formula
    audit_cols = [
        'study_id', 'first_author', 'year',
        'exposure_category', 'clock_type',
        'effect_type', 'effect_value', 'CI_lower', 'CI_upper', 'SE', 'SE_derived',
        'effect_class', 'flag_eaa_as_exposure', 'flag_wrong_outcome',
        'yi', 'sei', 'log_transformed',
        'poolable_primary', 'poolable_pool_B_std', 'poolable_pool_C_OR', 'poolable_pool_D_HR',
        'yi_sei_formula', 'conversion_notes',
    ]
    audit = out_final[[c for c in audit_cols if c in out_final.columns]].copy()
    audit.to_excel(writer, sheet_name='AUDIT_LOG', index=False)

print(f"\n✓ Saved: {OUTPUT_FILE}")
print("  Sheets: all_112_standardized | pool_A_primary_beta | pool_B_std_beta |")
print("          pool_C_OR_log | pool_D_HR_log | NEEDS_REVIEW_manual | AUDIT_LOG")
print("\nNext step: validate NEEDS_REVIEW sheet with professor, then run step7_meta_model.py")

import pandas as pd
import numpy as np

# ======================
# FILE INPUT
# ======================
INPUT_FILE = "extraction_raw_114studies.xlsx"

df = pd.read_excel(INPUT_FILE)

print("Rows loaded:", len(df))

# ======================
# NORMALIZZAZIONE STUDY DESIGN
# ======================
design_map = {
    "cross-sectional": "cross-sectional",
    "Cross-sectional": "cross-sectional",
    "cross-sectional observational": "cross-sectional",
    "cross-sectional observational study": "cross-sectional",
    "observational cross-sectional": "cross-sectional",

    "cohort study": "cohort",
    "prospective cohort": "cohort",
    "prospective cohort study": "cohort",
    "retrospective cohort study": "cohort",
    "longitudinal cohort": "cohort",
    "longitudinal cohort study": "cohort",
    "observational cohort study": "cohort",

    "case-control": "case-control",
    "nested case-control": "case-control",

    "Randomized controlled trial, 2x2x2 factorial design": "RCT",

    "Two-sample Mendelian randomization": "Mendelian randomization",
    "Two-sample bidirectional Mendelian randomization": "Mendelian randomization"
}

df["study_design_clean"] = df["study_design"].replace(design_map)

# ======================
# NORMALIZZAZIONE COUNTRY
# ======================
country_map = {
    "USA": "United States",
    "United States": "United States",
    "UK": "United Kingdom",
    "Scotland, UK": "United Kingdom"
}

df["country_clean"] = df["country"].replace(country_map)

# ======================
# DERIVAZIONE SE DA CI
# ======================
def derive_se(row):
    if pd.isna(row["SE"]) and not pd.isna(row["CI_lower"]) and not pd.isna(row["CI_upper"]):
        return (row["CI_upper"] - row["CI_lower"]) / 3.92
    return row["SE"]

df["SE_derived"] = df.apply(derive_se, axis=1)

# ======================
# REPORT MISSING
# ======================
missing_report = df[[
    "study_id",
    "first_author",
    "effect_value",
    "CI_lower",
    "CI_upper",
    "SE",
    "SE_derived",
    "p_value"
]].isna()

missing_summary = missing_report.sum()

missing_summary.to_excel("extraction_missing_report.xlsx")

# ======================
# SALVATAGGIO DATASET FINALE
# ======================
df.to_excel("extraction_raw_114studies_clean.xlsx", index=False)

print("Saved extraction_raw_114studies_clean.xlsx")
print("Saved extraction_missing_report.xlsx")

# ======================
# MINI TABLE 1
# ======================
print("\n==== STUDY DESIGN ====")
print(df["study_design_clean"].value_counts())

print("\n==== COUNTRY ====")
print(df["country_clean"].value_counts())

print("\n==== CLOCK TYPES ====")
print(df["clock_type"].value_counts().head(10))

print("\n==== SAMPLE SIZE ====")
print(df["sample_size"].describe())
import pandas as pd

# ===== FILE INPUT =====
INPUT_FILE = "step5_analytic_long_yi_sei.xlsx"

# ===== CARICAMENTO DATASET =====
df = pd.read_excel(INPUT_FILE)

print("Rows loaded:", len(df))

# ===== COLONNE DA TENERE (DATI GREZZI) =====
keep_cols = [
"study_id",
"first_author",
"year",
"journal",
"country",
"doi",
"study_design",
"sample_size",
"mean_age",
"exposure_category",
"exposure_variable",
"exposure_type",
"clock_type",
"EAA_metric",
"effect_type",
"effect_value",
"CI_lower",
"CI_upper",
"SE",
"p_value",
"adjusted_model",
"covariates"
]

df_raw = df[keep_cols].copy()

# ===== SALVATAGGIO TABELLA GREZZA =====
df_raw.to_excel("extraction_raw_114studies.xlsx", index=False)

print("Saved extraction_raw_114studies.xlsx")

# ===== TABLE 1 AUTOMATICA =====

print("\n===== SUMMARY TABLE =====")

print("\nStudies by design:")
print(df_raw["study_design"].value_counts())

print("\nStudies by clock:")
print(df_raw["clock_type"].value_counts())

print("\nExposure categories:")
print(df_raw["exposure_category"].value_counts())

print("\nCountry distribution:")
print(df_raw["country"].value_counts().head(10))

print("\nSample size statistics:")
print(df_raw["sample_size"].describe())

print("\nYear distribution:")
print(df_raw["year"].describe())
import pandas as pd
import re

INPUT_CSV = "csv-epigenetic-set.csv"          # metti il file nella stessa cartella dello script
OUTPUT_XLSX = "screening_step3.xlsx"

# Motivi standard (titolo/abstract) — userai questi per Exclude
REASONS = [
    "Wrong outcome (not epigenetic age/EAA)",
    "Wrong exposure (not modifiable determinant)",
    "Wrong population (not adults / animals / pediatric only)",
    "Wrong study type (review/editorial/commentary/protocol)",
    "Methodological only (clock development/validation without exposure-outcome)",
    "Other"
]

def norm_text(x: str) -> str:
    x = "" if pd.isna(x) else str(x)
    x = re.sub(r"\s+", " ", x).strip().lower()
    return x

df = pd.read_csv(INPUT_CSV)

# Normalizzazioni utili
df["title_norm"] = df["Title"].map(norm_text)
# Se non hai abstract nel CSV, lo screening si fa su titolo + (poi full-text/abstract da txt)
# Qui comunque prepariamo flags sul titolo.
title = df["title_norm"]

# Keyword flags (solo per velocizzare, NON per decidere automaticamente)
df["flag_clock"] = title.str.contains(r"\b(epigenetic|methylation|grim(age)?|phenoage|horvath|hannum|dunedin)\b", regex=True)
df["flag_acceleration"] = title.str.contains(r"\b(acceleration|age acceleration|eaa)\b", regex=True)
df["flag_modifiable"] = title.str.contains(r"\b(diet|nutrition|exercise|physical activity|smok|tobacco|alcohol|sleep|stress|pollution|air|socioeconomic|ses|occupation|lifestyle)\b", regex=True)
df["flag_review_like"] = title.str.contains(r"\b(review|systematic review|meta-analysis|protocol|editorial|commentary)\b", regex=True)

# Colonne screening (da compilare a mano)
df["decision"] = ""            # Include / Exclude / Unclear
df["reason_exclusion"] = ""    # scegli una voce da REASONS
df["notes"] = ""

# Ordine colonne (pulito)
cols_out = [
    "PMID", "DOI", "Publication Year", "Journal/Book", "First Author",
    "Title",
    "flag_clock", "flag_acceleration", "flag_modifiable", "flag_review_like",
    "decision", "reason_exclusion", "notes"
]
df_out = df[cols_out].sort_values(["flag_review_like", "flag_modifiable"], ascending=[True, False])

# Scrive Excel con 2 sheet: screening + reasons
with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
    df_out.to_excel(writer, index=False, sheet_name="screening")
    pd.DataFrame({"allowed_reasons": REASONS}).to_excel(writer, index=False, sheet_name="reasons")

print(f"Creato: {OUTPUT_XLSX}")
print("Ora apri l’Excel e compila 'decision' e (se Exclude) 'reason_exclusion'.")
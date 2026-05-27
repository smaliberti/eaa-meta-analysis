import pandas as pd
import re

# =========================
# CONFIG
# =========================
INPUT_CSV = "scopus_export.csv"   # <-- rinomina il tuo file così, oppure cambia qui il nome
OUTPUT_XLSX = "screening_step3_scopus.xlsx"

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

def first_author(authors: str) -> str:
    """
    Scopus 'Authors' spesso è tipo: "Rossi, M.; Bianchi, L.; ..."
    Prendiamo il primo blocco prima di ';'
    """
    if pd.isna(authors):
        return ""
    s = str(authors).strip()
    if ";" in s:
        return s.split(";", 1)[0].strip()
    return s

df = pd.read_csv(INPUT_CSV)

# Colonne attese Scopus:
# Authors | Title | Year | Source title | DOI | Abstract | Index Keywords | PubMed ID | EID
# (Author full names / Author(s) ID / Link possono esserci ma non sono obbligatorie)

# ID record robusto: usa PubMed ID se c'è, altrimenti EID
df["PubMed ID"] = df.get("PubMed ID", "")
df["EID"] = df.get("EID", "")
df["record_id"] = df["PubMed ID"].fillna("").astype(str).str.strip()
df.loc[df["record_id"].eq("") | df["record_id"].eq("nan"), "record_id"] = df["EID"].fillna("").astype(str).str.strip()

# Normalizzazioni
df["Title"] = df.get("Title", "")
df["title_norm"] = df["Title"].map(norm_text)

title = df["title_norm"]

# Flag keywords (solo per aiutare, NON decide in automatico)
df["flag_clock"] = title.str.contains(r"\b(epigenetic|methylation|dna methylation|grim(age)?|phenoage|horvath|hannum|dunedin|pace)\b", regex=True)
df["flag_acceleration"] = title.str.contains(r"\b(acceleration|age acceleration|eaa|pace of aging|dunedinpace)\b", regex=True)
df["flag_modifiable"] = title.str.contains(r"\b(diet|nutrition|exercise|physical activity|sedentary|smok|tobacco|alcohol|sleep|stress|pollution|air|pfas|phthalate|toxin|socioeconomic|ses|occupation|lifestyle)\b", regex=True)
df["flag_review_like"] = title.str.contains(r"\b(review|systematic review|meta-analysis|protocol|editorial|commentary|letter)\b", regex=True)

# Colonne decisione manuale
df["decision"] = ""            # Include / Exclude / Unclear
df["reason_exclusion"] = ""    # scegli da REASONS
df["notes"] = ""

# Campo autore e journal/anno coerenti con PubMed-style
df["First Author"] = df.get("Authors", "").map(first_author)
df["Publication Year"] = df.get("Year", "")
df["Journal/Book"] = df.get("Source title", "")

cols_out = [
    "record_id", "PubMed ID", "EID",
    "DOI", "Publication Year", "Journal/Book", "First Author",
    "Title",
    "flag_clock", "flag_acceleration", "flag_modifiable", "flag_review_like",
    "decision", "reason_exclusion", "notes"
]

# Mantieni solo le colonne che esistono davvero
cols_out = [c for c in cols_out if c in df.columns]

df_out = df[cols_out].sort_values(["flag_review_like", "flag_modifiable"], ascending=[True, False])

with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
    df_out.to_excel(writer, index=False, sheet_name="screening")
    pd.DataFrame({"allowed_reasons": REASONS}).to_excel(writer, index=False, sheet_name="reasons")

print(f"Creato: {OUTPUT_XLSX}")
print("Ora apri l’Excel e compila 'decision' e (se Exclude) 'reason_exclusion'.")




import pandas as pd
import re
from pathlib import Path

# =========================
# CONFIG
# =========================
CSV_PATH = Path("scopus_export.csv")  # <-- stesso nome del file usato nel code 1
SCREENING_XLSX_PATH = Path("screening_step3_scopus.xlsx")  # se esiste, integra decisioni
OUTPUT_XLSX_PATH = Path("screening_step3_scopus_autosuggest.xlsx")

# =========================
# HELPERS
# =========================
def norm(s: str) -> str:
    s = "" if s is None or pd.isna(s) else str(s)
    s = re.sub(r"\s+", " ", s).strip().lower()
    return s

def first_author(authors: str) -> str:
    if pd.isna(authors):
        return ""
    s = str(authors).strip()
    if ";" in s:
        return s.split(";", 1)[0].strip()
    return s

# =========================
# MAIN
# =========================
def main():
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Non trovo {CSV_PATH}. Mettilo nella cartella del progetto.")

    df = pd.read_csv(CSV_PATH)

    # Colonne base Scopus
    df["Title"] = df.get("Title", "")
    df["Abstract"] = df.get("Abstract", "")

    # ID robusto: PubMed ID -> altrimenti EID
    df["PubMed ID"] = df.get("PubMed ID", "")
    df["EID"] = df.get("EID", "")
    df["record_id"] = df["PubMed ID"].fillna("").astype(str).str.strip()
    df.loc[df["record_id"].eq("") | df["record_id"].eq("nan"), "record_id"] = df["EID"].fillna("").astype(str).str.strip()

    df["title_norm"] = df["Title"].map(norm)
    df["abstract_norm"] = df["Abstract"].map(norm)

    # Check: quanti abstract non vuoti?
    abs_len = df["abstract_norm"].str.len()
    print("Abstract length summary:")
    print(abs_len.describe())
    print(f"Abstract non-vuoti: {(abs_len > 0).sum()} / {len(df)}")

    # Testo combinato
    text = (df["title_norm"].fillna("") + " " + df["abstract_norm"].fillna("")).map(norm)

    # Pattern triage
    CLOCK_PAT = r"\b(?:epigenetic|dna methylation|methylation|epigenetic age|age acceleration|eaa|horvath|hannum|phenoage|grimage|dunedinpoam|dunedinpace|pace of aging)\b"
    MOD_PAT = r"\b(?:modifiable|lifestyle|behavioral|environmental|social|psychosocial|diet|nutrition|dietary|physical activity|exercise|sedentary|smok|tobacco|alcohol|sleep|stress|pollution|air pollution|pfas|phthalate|toxin|occupational|occupation|socioeconomic|ses)\b"
    REVIEW_PAT = r"\b(?:systematic review|meta-analysis|review|protocol|editorial|commentary|letter)\b"

    df["auto_has_clock"] = text.str.contains(CLOCK_PAT, regex=True)
    df["auto_has_modifiable"] = text.str.contains(MOD_PAT, regex=True)
    df["auto_review_like"] = df["title_norm"].str.contains(REVIEW_PAT, regex=True)

    def suggest(row):
        if row["auto_review_like"]:
            return ("Exclude", "Review/editorial/protocol (title)", "high")
        if not row["auto_has_clock"]:
            return ("Exclude", "No epigenetic clock/EAA terms in title+abstract", "high")
        if row["auto_has_clock"] and row["auto_has_modifiable"]:
            return ("Include", "Clock/EAA + modifiable determinant terms found", "high")
        if row["auto_has_clock"] and not row["auto_has_modifiable"]:
            return ("Unclear", "Clock/EAA found, modifiable determinant not obvious", "medium")
        return ("Unclear", "Insufficient information", "low")

    tmp = df.apply(suggest, axis=1, result_type="expand")
    df["auto_suggestion"] = tmp[0]
    df["auto_reason"] = tmp[1]
    df["auto_confidence"] = tmp[2]

    # Campi leggibili stile PubMed
    df["First Author"] = df.get("Authors", "").map(first_author)
    df["Publication Year"] = df.get("Year", "")
    df["Journal/Book"] = df.get("Source title", "")

    # Merge con screening manuale (se esiste)
    if SCREENING_XLSX_PATH.exists():
        prev = pd.read_excel(SCREENING_XLSX_PATH, sheet_name="screening")
        prev["record_id"] = prev.get("record_id", prev.get("PubMed ID", prev.get("EID", ""))).astype(str)

        keep_cols = ["record_id", "decision", "reason_exclusion", "notes"]
        for c in keep_cols:
            if c not in prev.columns:
                prev[c] = ""
        prev = prev[keep_cols]

        df = df.merge(prev, on="record_id", how="left")
    else:
        df["decision"] = ""
        df["reason_exclusion"] = ""
        df["notes"] = ""

    cols = [
        "record_id", "PubMed ID", "EID", "DOI",
        "Publication Year", "Journal/Book", "First Author",
        "Title",
        "auto_suggestion", "auto_confidence", "auto_reason",
        "decision", "reason_exclusion", "notes"
    ]
    cols = [c for c in cols if c in df.columns]
    df_out = df[cols].copy()

    # Ordine: Include -> Unclear -> Exclude
    order_map = {"Include": 0, "Unclear": 1, "Exclude": 2}
    df_out["auto_sort"] = df_out["auto_suggestion"].map(order_map).fillna(9).astype(int)
    df_out = df_out.sort_values(["auto_sort", "auto_confidence"]).drop(columns=["auto_sort"])

    with pd.ExcelWriter(OUTPUT_XLSX_PATH, engine="openpyxl") as w:
        df_out.to_excel(w, index=False, sheet_name="screening_autosuggest")

    print(f"\nSaved: {OUTPUT_XLSX_PATH}")
    print("Usa auto_suggestion come guida; compila decision manualmente (Include/Exclude/Unclear).")

if __name__ == "__main__":
    main()





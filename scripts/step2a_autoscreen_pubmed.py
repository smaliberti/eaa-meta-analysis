import pandas as pd
import re
from pathlib import Path

# =========================
# CONFIG
# =========================
CSV_PATH = Path("csv-epigenetic-set.csv")
ABSTRACT_TXT_PATH = Path("abstract-epigenetic-set.txt")
SCREENING_XLSX_PATH = Path("screening_step3.xlsx")  # se esiste, integra decisioni
OUTPUT_XLSX_PATH = Path("screening_step3_autosuggest.xlsx")


# =========================
# HELPERS
# =========================
def norm(s: str) -> str:
    s = "" if s is None else str(s)
    s = re.sub(r"\s+", " ", s).strip().lower()
    return s


def parse_abstracts_numbered_pubmed_export(txt: str) -> dict:
    """
    Parser per export PubMed in formato "numerato", tipo:

    1. Journal. date...
    Title.
    Authors.
    Author information:
    ...
    BACKGROUND: ... METHODS: ... RESULTS: ... CONCLUSIONS: ...
    ...
    DOI: ...
    PMCID: ...
    PMID: 12345678

    Restituisce: dict {PMID: abstract_norm}
    """
    # Split su inizio record: newline + numero + punto + spazio (lookahead)
    blocks = re.split(r"\n(?=\d+\.\s)", txt.strip())
    pmid_to_abs = {}

    for b in blocks:
        b = b.strip()
        if not b:
            continue

        # Estrai PMID (obbligatorio)
        m_pmid = re.search(r"\bPMID:\s*(\d+)", b)
        if not m_pmid:
            continue
        pmid = m_pmid.group(1)

        # Taglia via tutto dopo il PMID (pulizia)
        b_clean = re.split(r"\bPMID:\s*\d+\b", b, maxsplit=1)[0]

        # Caso A: abstract strutturato che inizia con BACKGROUND:
        m_struct = re.search(r"\b(BACKGROUND:.*)", b_clean, flags=re.DOTALL)
        if m_struct:
            abstract = m_struct.group(1)
        else:
            # Caso B: nessun BACKGROUND: -> prendi testo dopo "Author information:"
            if "Author information:" in b_clean:
                abstract = b_clean.split("Author information:", 1)[1]
            else:
                # Fallback: usa tutto (ultimo resort)
                abstract = b_clean

        # Tronca a DOI/PMCID se presenti
        abstract = re.split(r"\bDOI:\b|\bPMCID:\b", abstract, maxsplit=1)[0]

        pmid_to_abs[pmid] = norm(abstract)

    return pmid_to_abs


# =========================
# MAIN
# =========================
def main():
    # 1) Load metadata CSV
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Non trovo {CSV_PATH}. Mettilo nella cartella del progetto.")

    df = pd.read_csv(CSV_PATH)
    if "PMID" not in df.columns or "Title" not in df.columns:
        raise ValueError("Il CSV deve contenere almeno le colonne 'PMID' e 'Title'.")

    df["PMID"] = df["PMID"].astype(str)
    df["title_norm"] = df["Title"].map(norm)

    # 2) Parse abstracts from txt (formato numerato)
    if not ABSTRACT_TXT_PATH.exists():
        raise FileNotFoundError(f"Non trovo {ABSTRACT_TXT_PATH}. Mettilo nella cartella del progetto.")

    txt = ABSTRACT_TXT_PATH.read_text(encoding="utf-8", errors="ignore")
    pmid_to_abs = parse_abstracts_numbered_pubmed_export(txt)

    df["abstract_norm"] = df["PMID"].map(lambda x: pmid_to_abs.get(str(x), ""))

    # 3) Verifica rapida: quanti abstract abbiamo agganciato?
    abs_len = df["abstract_norm"].str.len()
    print("Abstract length summary:")
    print(abs_len.describe())
    print(f"Abstract non-vuoti: {(abs_len > 0).sum()} / {len(df)}")

    # 4) Combined text for keyword triage
    text = (df["title_norm"].fillna("") + " " + df["abstract_norm"].fillna("")).map(norm)

    # 5) Keyword patterns (non-capturing)
    CLOCK_PAT = r"\b(?:epigenetic|methylation|dna methylation|epigenetic age|eaa|age acceleration|horvath|hannum|phenoage|grimage|dunedinpoam|dunedinpace)\b"
    MOD_PAT = r"\b(?:modifiable|lifestyle|behavioral|environmental|social|psychosocial|diet|nutrition|dietary pattern|physical activity|exercise|sedentary|smok|tobacco|alcohol|sleep|stress|pollution|air pollution|toxin|socioeconomic|ses|occupation|occupational)\b"
    REVIEW_PAT = r"\b(?:systematic review|meta-analysis|review|protocol|editorial|commentary|letter)\b"

    df["auto_has_clock"] = text.str.contains(CLOCK_PAT, regex=True)
    df["auto_has_modifiable"] = text.str.contains(MOD_PAT, regex=True)
    df["auto_review_like"] = df["title_norm"].str.contains(REVIEW_PAT, regex=True)

    # 6) Suggestion rules (triage, not final decision)
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

    # 7) Merge your manual screening columns if present
    if SCREENING_XLSX_PATH.exists():
        prev = pd.read_excel(SCREENING_XLSX_PATH, sheet_name="screening")
        prev["PMID"] = prev["PMID"].astype(str)

        keep_cols = ["PMID", "decision", "reason_exclusion", "notes"]
        for c in keep_cols:
            if c not in prev.columns:
                prev[c] = ""

        prev = prev[keep_cols]
        df = df.merge(prev, on="PMID", how="left")
    else:
        df["decision"] = ""
        df["reason_exclusion"] = ""
        df["notes"] = ""

    # 8) Output Excel
    # Colonne "presentabili"
    cols = [
        "PMID", "DOI", "Publication Year", "Journal/Book", "First Author", "Title",
        "auto_suggestion", "auto_confidence", "auto_reason",
        "decision", "reason_exclusion", "notes"
    ]
    # Alcuni CSV potrebbero non avere tutte queste colonne: gestiamo in modo robusto
    existing_cols = [c for c in cols if c in df.columns]
    df_out = df[existing_cols].copy()

    # Ordina per priorità: Include -> Unclear -> Exclude
    order_map = {"Include": 0, "Unclear": 1, "Exclude": 2}
    df_out["auto_sort"] = df["auto_suggestion"].map(order_map).fillna(9).astype(int)
    df_out = df_out.sort_values(["auto_sort", "auto_confidence"]).drop(columns=["auto_sort"])

    with pd.ExcelWriter(OUTPUT_XLSX_PATH, engine="openpyxl") as w:
        df_out.to_excel(w, index=False, sheet_name="screening_autosuggest")

    print(f"\nSaved: {OUTPUT_XLSX_PATH}")
    print("Usa auto_suggestion come guida; compila decision manualmente (Include/Exclude/Unclear).")


if __name__ == "__main__":
    main()




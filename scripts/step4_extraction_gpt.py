"""
PDF -> OpenAI -> JSON -> append to extraction_master.xlsx

Prerequisiti:
pip install openai pymupdf pandas openpyxl

Imposta OPENAI_API_KEY come variabile ambiente (consigliato).
"""

import os
import re
import json
import time
from pathlib import Path

import fitz  # PyMuPDF
import pandas as pd
from openai import OpenAI


# =========================
# CONFIG
# =========================
MODEL = "gpt-4.1-mini"   # se vuoi cambiare, dimmelo e ti do alternative
MAX_PAGES = 10            # quante pagine leggere (puoi alzare a 10-12 se serve)
MAX_CHARS = 35000        # taglio di sicurezza del testo per evitare prompt enormi
SLEEP_BETWEEN_CALLS = 0.7

# >>>>>>> MODIFICA QUI I PERCORSI <<<<<<<
# ESEMPIO: singolo PDF
PDF_PATH = Path("/Users/pier/Desktop/pubUNISALERNO/Meta_EAA_Project/01_PDF_fulltext/Scopus/The-association-between-phenotypic-age-acceleration-and-the-risk-of-allcause-and-cancer-mortality-among-cancer-survivors-NHANES-19992018_2026_Nature-Research.pdf")

# Oppure: cartella di PDF (batch)
PDF_FOLDER = Path("/Users/pier/Desktop/pubUNISALERNO/Meta_EAA_Project/01_PDF_fulltext/Scopus")

# Il tuo master excel
MASTER_PATH = Path("/Users/pier/Desktop/pubUNISALERNO/Meta_EAA_Project/02_Data_extraction/extraction_master.xlsx")

# Se vuoi batch: True = scansiona tutta la cartella, False = usa PDF_PATH
RUN_BATCH = False


# =========================
# HELPERS
# =========================
def extract_text(pdf_path: Path, max_pages: int = MAX_PAGES) -> str:
    doc = fitz.open(pdf_path)
    chunks = []
    for i, page in enumerate(doc):
        if i >= max_pages:
            break
        chunks.append(page.get_text("text"))
    doc.close()
    text = "\n".join(chunks)

    if len(text) > MAX_CHARS:
        text = text[:MAX_CHARS]
    return text


def parse_json_loose(s: str) -> dict:
    """
    Gestisce output tipo:
    ```json
    {...}
    ```
    oppure testo extra prima/dopo.
    """
    if not s:
        raise ValueError("Empty model output")

    s = s.strip()

    # rimuovi fence ```json ... ``` / ``` ... ```
    s = re.sub(r"^```(?:json)?\s*", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s*```$", "", s)

    # se non parte con { prova a estrarre il primo {...}
    if not s.lstrip().startswith("{"):
        m = re.search(r"\{.*\}", s, flags=re.DOTALL)
        if not m:
            raise ValueError(f"Could not find JSON object in output. Head:\n{s[:500]}")
        s = m.group(0)

    return json.loads(s)


def build_prompt(paper_text: str, file_name: str, columns: list[str]) -> str:
    fields = "\n".join([f"- {c}" for c in columns])
    return f"""
You are extracting structured data for a meta-analysis on modifiable determinants and epigenetic age acceleration.

Return ONLY a JSON object with EXACTLY these keys (no extra keys):
{fields}

Rules:
- If a field is not reported, use null.
- Be conservative: do NOT guess.
- If multiple exposures/clocks exist, summarize in a compact, readable way (e.g., "MVPA; steps/day; sedentary time") rather than long lists.
- effect_value should be a SINGLE main estimate relevant for meta-analysis if possible; otherwise set null and describe in notes-like fields (covariates/adjusted_model/exposure_variable) what is available.
- p_value: numeric if possible; if only thresholds reported (e.g., "<0.001"), use 0.001.

The PDF file name is: {file_name}

Paper text (partial):
{paper_text}
""".strip()


def call_openai_json(client: OpenAI, prompt: str) -> dict:
    """
    Usa chat.completions per massima compatibilità.
    """
    # “JSON-mode” compatibile: chiediamo JSON esplicito nel prompt + parser loose
    resp = client.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "system", "content": "You extract information from scientific papers and output valid JSON only."},
            {"role": "user", "content": prompt},
        ],
        temperature=0,
    )

    out = resp.choices[0].message.content
    return parse_json_loose(out)


def ensure_master_exists(master_path: Path) -> pd.DataFrame:
    if master_path.exists():
        return pd.read_excel(master_path)

    # se non esiste, lo creiamo con le intestazioni standard
    cols = [
        "study_id","first_author","year","journal","country","doi","study_design","sample_size",
        "mean_age","percent_female","exposure_category","exposure_variable","exposure_type",
        "clock_type","EAA_metric","effect_type","effect_value","CI_lower","CI_upper","SE","p_value",
        "adjusted_model","covariates","conversion_method","yi","sei","direction_checked"
    ]
    df = pd.DataFrame(columns=cols)
    master_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(master_path, index=False)
    return df


def append_row(master_path: Path, row: dict) -> None:
    df_master = ensure_master_exists(master_path)
    cols = list(df_master.columns)

    df_new = pd.DataFrame([row]).reindex(columns=cols)
    df_out = pd.concat([df_master, df_new], ignore_index=True)
    df_out.to_excel(master_path, index=False)


# =========================
# MAIN
# =========================
def run_one(pdf_path: Path, client: OpenAI):
    df_master = ensure_master_exists(MASTER_PATH)
    columns = list(df_master.columns)

    text = extract_text(pdf_path, MAX_PAGES)
    prompt = build_prompt(text, pdf_path.name, columns)

    data = call_openai_json(client, prompt)

    # safety: assicurati che le chiavi corrispondano al master
    missing = [c for c in columns if c not in data]
    extra = [k for k in data.keys() if k not in columns]
    if missing:
        # metti null per ciò che manca
        for m in missing:
            data[m] = None
    if extra:
        # elimina chiavi extra
        for e in extra:
            data.pop(e, None)

    append_row(MASTER_PATH, data)
    print(f"✅ Appended: {pdf_path.name}")


def main():
    # API KEY da env
    if not os.getenv("OPENAI_API_KEY"):
        raise RuntimeError(
            "OPENAI_API_KEY non trovata. Impostala come variabile ambiente "
            "(PyCharm: Run/Debug Configurations -> Environment variables)."
        )

    client = OpenAI()

    if RUN_BATCH:
        pdfs = sorted(PDF_FOLDER.glob("*.pdf"))
        if not pdfs:
            raise FileNotFoundError(f"Nessun PDF trovato in: {PDF_FOLDER}")
        print(f"Trovati {len(pdfs)} PDF in {PDF_FOLDER}")

        for i, pdf in enumerate(pdfs, start=1):
            try:
                print(f"\n[{i}/{len(pdfs)}] Processing: {pdf.name}")
                run_one(pdf, client)
                time.sleep(SLEEP_BETWEEN_CALLS)
            except Exception as e:
                print(f"❌ ERRORE su {pdf.name}: {e}")
                # continua con i successivi
                continue
    else:
        if not PDF_PATH.exists():
            raise FileNotFoundError(f"PDF non trovato: {PDF_PATH}")
        run_one(PDF_PATH, client)


if __name__ == "__main__":
    main()
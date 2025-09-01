#!/usr/bin/env python3
"""
Gene/Protein mention counter for multiple PDFs — **No OCR**
Offline SciSpaCy NER + optional HGNC synonym normalization + optional whitelist.

Quick start (in your venv):
    # 1) Core deps
    python -m pip install --upgrade pip
    python -m pip install pymupdf tqdm "spacy==3.7.5"

    # 2) SciSpaCy without heavy extras (avoids nmslib build on Windows)
    python -m pip install "scispacy==0.5.5" --no-deps
    python -m pip install numpy==1.26.4 scikit-learn==1.7.1 joblib==1.5.1 pysbd==0.3.4 conllu==6.0.0

    # 3) SciSpaCy model wheel (JNLPBA)
    python -m pip install \
      https://github.com/allenai/scispacy/releases/download/v0.5.5/en_ner_jnlpba_md-0.5.5-py3-none-any.whl

Usage:
    python gene_counter_no_ocr.py "Papers/*.pdf" \
        --hgnc hgnc_complete_set.txt \
        --whitelist ecoli_symbols.txt \
        --global-out global_counts.csv \
        --per-file-out per_file_counts.csv \
        --top 30 \
        --n-process 1 \
        --ref-weight 0.5

Notes:
- No OCR: scanned PDFs (no text layer) will be skipped with a warning.
- Labels used: {"gene", "protein", "gene_or_gene_product"}.
- If you don't have HGNC or a species whitelist, just omit those flags.
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Tuple

import fitz  # PyMuPDF
from tqdm import tqdm

# Lazy import with a helpful message
try:
    import spacy
except ImportError:
    sys.exit("[ERROR] SciSpaCy/spaCy not installed. See Quick start at top of file.")

MODEL_NAME = "en_ner_jnlpba_md"  # SciSpaCy NER model (installed via wheel above)

# --------------------------------------------------------------------------------------
# PDF → text (no OCR)
# --------------------------------------------------------------------------------------

def extract_pdf_pages(pdf_path: pathlib.Path) -> List[str]:
    """Return a list of plain-text strings, one per page. Empty pages become ''."""
    pages: List[str] = []
    try:
        with fitz.open(pdf_path) as doc:
            for page in doc:
                txt = page.get_text("text") or ""
                pages.append(txt)
    except Exception as e:
        raise RuntimeError(f"{pdf_path.name}: {e}")
    return pages


def normalize_page_text(text: str) -> str:
    """Light normalization for NER friendliness: join hyphenated line-breaks, collapse spaces."""
    # Fix hyphenation across line breaks: e.g., "gene-\nname" -> "genename"
    text = re.sub(r"(\w)-\n(\w)", r"\1\2", text)
    # Remove stray line breaks while keeping paragraph separation (very light touch)
    text = re.sub(r"\n{2,}", "\n\n", text)
    text = text.replace("\n", " ")
    # Normalize multiple spaces
    text = re.sub(r"\s{2,}", " ", text).strip()
    return text


# --------------------------------------------------------------------------------------
# Lexicon helpers (HGNC + optional whitelist)
# --------------------------------------------------------------------------------------

def load_hgnc(tsv: pathlib.Path | None) -> Dict[str, str]:
    """Load HGNC TSV and build UPPERCASE alias→canonical map. Safe to call with None."""
    mapping: Dict[str, str] = {}
    if not tsv:
        return mapping
    if not tsv.exists():
        print(f"[WARN] HGNC file not found: {tsv}")
        return mapping

    import csv as _csv
    with tsv.open(encoding="utf-8") as fh:
        reader = _csv.DictReader(fh, delimiter="\t")
        required = {"symbol", "alias_symbol", "prev_symbol"}
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            print(f"[WARN] HGNC TSV missing columns: {missing}")
        for row in reader:
            symbol = (row.get("symbol") or "").upper()
            if not symbol:
                continue
            mapping[symbol] = symbol
            for col in ("alias_symbol", "prev_symbol"):
                val = row.get(col) or ""
                for token in (val.split("|") if val else []):
                    token = token.strip().upper()
                    if token:
                        mapping[token] = symbol
    return mapping


def load_whitelist(path: pathlib.Path | None) -> set[str]:
    """Load a simple newline-delimited whitelist of valid gene symbols (UPPERCASE)."""
    if not path:
        return set()
    if not path.exists():
        print(f"[WARN] whitelist not found: {path}")
        return set()
    syms = set()
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            syms.add(s.upper())
    return syms


# --------------------------------------------------------------------------------------
# Counting
# --------------------------------------------------------------------------------------

VALID_LABELS = {"gene", "protein", "gene_or_gene_product"}


def count_entities_in_pages(
    pages: List[str],
    nlp,
    hgnc_map: Dict[str, str],
    whitelist: set[str],
    ref_weight: float = 0.5,
    n_process: int = 1,
    batch_size: int = 1000,
) -> Counter:
    """Run NER per page with optional down-weighting for pages after 'References'."""
    cleaned: List[str] = [normalize_page_text(p) for p in pages]

    # Determine per-page weights: after first occurrence of 'References' heading → ref_weight
    weights: List[float] = []
    seen_refs = False
    for p in pages:
        if not seen_refs and re.search(r"\bReferences\b", p, flags=re.IGNORECASE):
            seen_refs = True
        weights.append(ref_weight if seen_refs else 1.0)

    counts: Counter = Counter()
    # Use nlp.pipe for efficiency
    for doc, w in zip(
        nlp.pipe(cleaned, n_process=max(1, n_process), batch_size=batch_size, disable=["parser", "attribute_ruler", "lemmatizer"]),
        weights,
    ):
        for ent in doc.ents:
            if ent.label_.lower() in VALID_LABELS:
                sym = ent.text.strip().upper()
                # Basic cleanup: collapse spaces and punctuation common in entities
                sym = re.sub(r"\s+", " ", sym)
                sym = sym.strip(";,:.()[]{}")
                # HGNC normalization first (if provided)
                sym = hgnc_map.get(sym, sym)
                # Whitelist filter (if provided)
                if whitelist and sym not in whitelist:
                    continue
                # Count with weight (Counter supports ints; we store scaled ints by 1000 to keep CSV neat)
                counts[sym] += w
    return counts


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Count Gene/Protein mentions across PDFs (No OCR)")
    ap.add_argument("pdfs", nargs="+", help="PDF glob(s), e.g. 'Papers/*.pdf' or 'A/*.pdf' 'B/*.pdf'")
    ap.add_argument("--hgnc", type=pathlib.Path, default=None, help="HGNC TSV for synonym merge (optional)")
    ap.add_argument("--whitelist", type=pathlib.Path, default=None, help="UPPERCASE gene symbols (one per line)")
    ap.add_argument("--global-out", type=pathlib.Path, default=pathlib.Path("global_counts.csv"))
    ap.add_argument("--per-file-out", type=pathlib.Path, default=pathlib.Path("per_file_counts.csv"))
    ap.add_argument("--model", default=MODEL_NAME, help=f"spaCy/SciSpaCy model (default: {MODEL_NAME})")
    ap.add_argument("--top", type=int, default=0, help="Print top-N genes to stdout")
    ap.add_argument("--n-process", type=int, default=1, help="spaCy n_process (>=1)")
    ap.add_argument("--batch-size", type=int, default=1000, help="spaCy pipe batch size")
    ap.add_argument("--ref-weight", type=float, default=0.5, help="Weight for pages after 'References' heading")
    return ap.parse_args()


def resolve_globs(patterns: Iterable[str]) -> List[pathlib.Path]:
    files: List[pathlib.Path] = []
    for patt in patterns:
        p = pathlib.Path(patt)
        if p.exists():             # 파일/폴더 경로를 그대로 받은 경우
            if p.is_file():
                files.append(p)
            else:
                files.extend(p.rglob("*.pdf"))
            continue
        # 글롭 패턴인 경우 (예: Papers/*.pdf)
        files.extend(pathlib.Path().glob(patt))
    # 중복 제거 (순서 유지)
    seen = set()
    unique: List[pathlib.Path] = []
    for f in files:
        r = f.resolve()
        if r not in seen:
            seen.add(r)
            unique.append(r)
    return unique


def write_global_counts(global_counts: Counter, out_path: pathlib.Path) -> None:
    out_path.write_text("gene,count\n", encoding="utf-8")
    with out_path.open("a", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        for gene, cnt in global_counts.most_common():
            # cnt may be float if ref-weight != 1; format with up to 3 decimals
            if isinstance(cnt, float) and not cnt.is_integer():
                w.writerow([gene, f"{cnt:.3f}"])
            else:
                w.writerow([gene, int(cnt)])


def write_per_file(per_file_counts: Dict[str, Counter], out_path: pathlib.Path) -> None:
    out_path.write_text("pdf,gene,count\n", encoding="utf-8")
    with out_path.open("a", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        for pdf_name, counter in per_file_counts.items():
            for gene, cnt in counter.most_common():
                if isinstance(cnt, float) and not cnt.is_integer():
                    w.writerow([pdf_name, gene, f"{cnt:.3f}"])
                else:
                    w.writerow([pdf_name, gene, int(cnt)])


def main() -> None:
    args = parse_args()

    # Load model
    try:
        nlp = spacy.load(args.model)
    except Exception as e:
        sys.exit(
            "[ERROR] Failed to load model. Ensure you installed the SciSpaCy wheel shown in the header.\n"
            f"        Model: {args.model}\n        Reason: {e}"
        )

    hgnc_map = load_hgnc(args.hgnc)
    if hgnc_map:
        print(f"[INFO] HGNC synonyms loaded: {len(hgnc_map):,} entries")
    whitelist = load_whitelist(args.whitelist)
    if not whitelist and hgnc_map:
        whitelist = set(hgnc_map.values())
        print(f"[INFO] Whitelist(auto): using {len(whitelist):,} HGNC symbols (keep only valid gene symbols)")
    if whitelist:
        print(f"[INFO] Whitelist loaded: {len(whitelist):,} symbols")

    pdf_files = resolve_globs(args.pdfs)
    if not pdf_files:
        sys.exit("[ERROR] No PDF files matched the given pattern(s).")

    global_counts: Counter = Counter()
    per_file: Dict[str, Counter] = {}

    for pdf in tqdm(pdf_files, desc="PDFs → NER"):
        try:
            pages = extract_pdf_pages(pdf)
            if not any(p.strip() for p in pages):
                print(f"[WARN] No text layer found (likely scanned): {pdf.name} — skipped (No OCR mode)")
                continue
            c = count_entities_in_pages(
                pages,
                nlp,
                hgnc_map=hgnc_map,
                whitelist=whitelist,
                ref_weight=args.ref_weight,
                n_process=max(1, args.n_process),
                batch_size=max(1, args.batch_size),
            )
            per_file[pdf.name] = c
            global_counts.update(c)
        except Exception as e:
            print(f"[WARN] {pdf.name}: {e}")

    write_global_counts(global_counts, args.global_out)
    write_per_file(per_file, args.per_file_out)

    print(f"[INFO] Saved global counts → {args.global_out.resolve()}")
    print(f"[INFO] Saved per-file counts → {args.per_file_out.resolve()}")

    if args.top > 0:
        print("\nTop genes (global):")
        for gene, cnt in global_counts.most_common(args.top):
            cnt_disp = f"{cnt:.3f}" if isinstance(cnt, float) and not cnt.is_integer() else str(int(cnt))
            print(f"{gene:>12} : {cnt_disp}")


if __name__ == "__main__":
    main()

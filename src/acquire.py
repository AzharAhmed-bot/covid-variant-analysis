"""Fetch SARS-CoV-2 spike protein sequences with variant, place and date.

The project needs, for each isolate: an accession, a Pango lineage, a country
and geographic region, a collection date, and the spike protein sequence.

Two routes are supported.

  api   NCBI Datasets v2 virus API. Fully programmatic, so the notebook runs
        top to bottom without manual steps. Requires outbound access to
        api.ncbi.nlm.nih.gov.

  local A metadata table plus a protein FASTA exported from the NCBI Virus
        portal. Used when the API is unreachable.

Both routes produce the same canonical table, defined by CANONICAL_COLUMNS.
"""
from __future__ import annotations

import io
import json
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
SARS2_TAXON = "2697049"
USER_AGENT = "strathmore-bioinformatics-coursework/1.0"

CANONICAL_COLUMNS = [
    "accession",        # str, NCBI nucleotide accession, unique
    "lineage",          # str, Pango lineage as reported by NCBI
    "variant",          # str, Alpha / Delta / Omicron (derived from lineage)
    "country",          # str
    "region",           # str, continent-level geographic region
    "collection_date",  # datetime64, isolate collection date
    "spike",            # str, spike protein sequence as deposited
]


def _request(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                               "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def api_dataset_report(
    lineage: str,
    page_size: int = 100,
    max_records: int = 500,
    released_since: str | None = None,
    complete_only: bool = True,
    throttle: float = 0.4,
) -> list[dict]:
    """Page through the virus dataset report for one Pango lineage.

    Returns the raw report dictionaries so that field names can be inspected
    before mapping. Kept deliberately thin: NCBI has changed field names
    between releases, so mapping happens in `report_to_frame`.
    """
    records: list[dict] = []
    page_token = None
    while len(records) < max_records:
        params = {
            "filter.pangolin_classification": lineage,
            "page_size": str(min(page_size, max_records - len(records))),
        }
        if complete_only:
            params["filter.complete_only"] = "true"
        if released_since:
            params["filter.released_since"] = released_since
        if page_token:
            params["page_token"] = page_token
        url = (f"{DATASETS_BASE}/virus/taxon/{SARS2_TAXON}/genome/dataset_report"
               f"?{urllib.parse.urlencode(params)}")
        payload = json.loads(_request(url))
        batch = payload.get("reports") or []
        if not batch:
            break
        records.extend(batch)
        page_token = payload.get("next_page_token")
        if not page_token:
            break
        time.sleep(throttle)
    return records[:max_records]


def api_download_proteins(accessions: list[str], out_zip: str | Path,
                          chunk: int = 200, throttle: float = 0.5) -> Path:
    """Download annotated protein FASTA for the given accessions."""
    out_zip = Path(out_zip)
    out_zip.parent.mkdir(parents=True, exist_ok=True)
    parts = []
    for i in range(0, len(accessions), chunk):
        batch = accessions[i : i + chunk]
        url = (f"{DATASETS_BASE}/virus/accession/{','.join(batch)}/download"
               f"?include_annotation_type=PROT_FASTA")
        parts.append(_request(url, timeout=600))
        time.sleep(throttle)
    # Concatenating zips is not valid; keep them separate on disk instead.
    if len(parts) == 1:
        out_zip.write_bytes(parts[0])
    else:
        with zipfile.ZipFile(out_zip, "w") as combined:
            for n, raw in enumerate(parts):
                with zipfile.ZipFile(io.BytesIO(raw)) as z:
                    for name in z.namelist():
                        if name.endswith(".faa") or name.endswith(".fasta"):
                            combined.writestr(f"part{n}_{Path(name).name}",
                                              z.read(name))
    return out_zip


def extract_spike_from_protein_fasta(raw: bytes) -> dict[str, str]:
    """Pull surface glycoprotein records out of a protein FASTA zip.

    Returns {accession: spike sequence}. Records are matched on the product
    name in the defline, accepting the several spellings NCBI has used.
    """
    spike_names = ("surface glycoprotein", "spike glycoprotein", "spike protein")
    out: dict[str, str] = {}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for name in z.namelist():
            if not (name.endswith(".faa") or name.endswith(".fasta")):
                continue
            text = z.read(name).decode("utf-8", errors="replace")
            header, chunks = None, []
            for line in text.splitlines():
                if line.startswith(">"):
                    if header is not None:
                        _store(out, header, "".join(chunks), spike_names)
                    header, chunks = line[1:], []
                else:
                    chunks.append(line.strip())
            if header is not None:
                _store(out, header, "".join(chunks), spike_names)
    return out


def _store(out: dict, header: str, seq: str, spike_names: tuple[str, ...]) -> None:
    low = header.lower()
    if not any(n in low for n in spike_names):
        return
    # Deflines carry the source nucleotide accession in brackets or as the
    # first whitespace-delimited token; try both.
    acc = None
    if "[" in header and "]" in header:
        for field in header.split("["):
            val = field.split("]")[0].strip()
            if val.count(".") == 1 and len(val) < 20 and val[0].isalpha():
                acc = val
    if acc is None:
        acc = header.split()[0].split("|")[-1]
    if acc and seq:
        out[acc] = seq.upper()


def load_local(metadata_path: str | Path, fasta_path: str | Path) -> pd.DataFrame:
    """Load an NCBI Virus portal export: metadata table plus protein FASTA."""
    meta = pd.read_csv(metadata_path)
    seqs: dict[str, str] = {}
    header, chunks = None, []
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(chunks).upper()
                header, chunks = line[1:].split()[0].split("|")[0], []
            else:
                chunks.append(line.strip())
    if header:
        seqs[header] = "".join(chunks).upper()
    meta["spike"] = meta.iloc[:, 0].astype(str).map(seqs)
    return meta

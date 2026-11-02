"""Deepen the sample for chosen lineages and regions."""
from __future__ import annotations

import collections
import csv
import json
import sys
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import build_dataset as bd
import spike_reference as ref

WANTED_REGIONS = {"Europe", "North America", "Africa", "Asia", "South America"}
METADATA_PATH = Path("data/raw/metadata.csv")
FASTA_PATH = Path("data/raw/spike_proteins.fasta")


def collect_new_records(lineage: str, cap_per_region: int, max_pages: int,
                        already_held: set[str]) -> list[dict]:
    """Page the dataset report for one lineage, keeping unseen records."""
    kept: list[dict] = []
    per_region: collections.Counter = collections.Counter()
    token, scanned = None, 0

    for page in range(max_pages):
        params = {"filter.pangolin_classification": lineage,
                  "page_size": bd.PAGE_SIZE}
        if token:
            params["page_token"] = token
        try:
            payload = json.loads(bd.fetch(
                f"{bd.DATASETS}/virus/taxon/sars2/dataset_report?"
                + urllib.parse.urlencode(params)))
        except Exception as exc:
            print(f"  {lineage} page {page} failed: {exc}", flush=True)
            break

        reports = payload.get("reports") or []
        if not reports:
            break

        for rec in reports:
            scanned += 1
            iso = rec.get("isolate") or {}
            loc = rec.get("location") or {}
            region = loc.get("geographic_region")
            date = iso.get("collection_date")
            if (not region or not date or region not in WANTED_REGIONS
                    or rec.get("is_annotated") is False
                    or rec["accession"] in already_held
                    or per_region[region] >= cap_per_region):
                continue
            kept.append({
                "accession": rec["accession"],
                "lineage": lineage,
                "variant": ref.variant_from_pango(lineage),
                "country": (loc.get("geographic_location") or "").split(":")[0].strip(),
                "region": region,
                "collection_date": date,
                "genome_length": rec.get("length"),
                "completeness": rec.get("completeness"),
            })
            per_region[region] += 1
            already_held.add(rec["accession"])

        token = payload.get("next_page_token")
        if page % 20 == 0:
            print(f"  {lineage} page {page:3d} scanned {scanned:6d} "
                  f"new {len(kept):5d} {dict(per_region)}", flush=True)
        if not token:
            break
        time.sleep(bd.DATASETS_PAUSE)

    print(f"{lineage}: scanned {scanned}, new {len(kept)} {dict(per_region)}",
          flush=True)
    return kept


def append_records(rows: list[dict], sequences: dict[str, str]) -> None:
    """Append new records to the metadata table and the sequence file."""
    with open(METADATA_PATH, "a", newline="") as handle:
        csv.DictWriter(handle, fieldnames=list(rows[0].keys())).writerows(rows)
    with open(FASTA_PATH, "a") as handle:
        for row in rows:
            handle.write(f">{row['accession']}|{row['lineage']}|{row['variant']}|"
                         f"{row['region']}|{row['collection_date']}\n")
            seq = sequences[row["accession"]]
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + "\n")


def main(cap_per_region: int, max_pages: int, lineages: list[str]) -> None:
    already_held = {r["accession"] for r in csv.DictReader(open(METADATA_PATH))}
    print(f"{len(already_held)} accessions already held", flush=True)

    for lineage in lineages:
        kept = collect_new_records(lineage, cap_per_region, max_pages, already_held)
        if not kept:
            continue
        sequences = bd.fetch_spikes([r["accession"] for r in kept])
        rows = [r for r in kept if r["accession"] in sequences]
        print(f"{lineage}: {len(rows)} with spike", flush=True)
        if rows:
            append_records(rows, sequences)
            print(f"{lineage}: appended", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit("usage: topup.py <cap_per_region> <max_pages> <lineage> [lineage ...]")
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3:])

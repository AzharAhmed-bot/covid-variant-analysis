"""Deepen the sample for chosen lineages and regions.

Within a variant, spike carries limited diversity: thousands of isolates
collapse to a few hundred distinct mutation profiles. The geographic task is
therefore limited by distinct profiles rather than by record count, and the
only way to gain power is to scan more records and keep the new profiles.

Run: python3 -I src/topup.py <cap_per_region> <max_pages> <lineage> [lineage ...]
"""
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

CAP = int(sys.argv[1])
PAGES = int(sys.argv[2])
LINEAGES = sys.argv[3:]
# Regions worth deepening: the ones with enough records to support a model.
WANTED = {"Europe", "North America", "Africa", "Asia", "South America"}

meta_path = Path("data/raw/metadata.csv")
fasta_path = Path("data/raw/spike_proteins.fasta")
existing = {r["accession"] for r in csv.DictReader(open(meta_path))}
print(f"{len(existing)} accessions already held", flush=True)

for lineage in LINEAGES:
    kept: list[dict] = []
    per_region: collections.Counter = collections.Counter()
    token, scanned = None, 0
    for page in range(PAGES):
        params = {"filter.pangolin_classification": lineage, "page_size": bd.PAGE_SIZE}
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
            iso, loc = rec.get("isolate") or {}, rec.get("location") or {}
            region, date = loc.get("geographic_region"), iso.get("collection_date")
            if (not region or not date or region not in WANTED
                    or rec.get("is_annotated") is False
                    or rec["accession"] in existing
                    or per_region[region] >= CAP):
                continue
            kept.append({
                "accession": rec["accession"], "lineage": lineage,
                "variant": ref.variant_from_pango(lineage),
                "country": (loc.get("geographic_location") or "").split(":")[0].strip(),
                "region": region, "collection_date": date,
                "genome_length": rec.get("length"),
                "completeness": rec.get("completeness"),
            })
            per_region[region] += 1
            existing.add(rec["accession"])
        token = payload.get("next_page_token")
        if page % 20 == 0:
            print(f"  {lineage} page {page:3d} scanned {scanned:6d} "
                  f"new {len(kept):5d} {dict(per_region)}", flush=True)
        if not token:
            break
        time.sleep(bd.DATASETS_PAUSE)

    print(f"{lineage}: scanned {scanned}, new {len(kept)} {dict(per_region)}", flush=True)
    if not kept:
        continue
    seqs = bd.fetch_spikes([r["accession"] for r in kept])
    rows = [r for r in kept if r["accession"] in seqs]
    print(f"{lineage}: {len(rows)} with spike", flush=True)
    if not rows:
        continue
    with open(meta_path, "a", newline="") as handle:
        csv.DictWriter(handle, fieldnames=list(rows[0].keys())).writerows(rows)
    with open(fasta_path, "a") as handle:
        for r in rows:
            handle.write(f">{r['accession']}|{r['lineage']}|{r['variant']}|"
                         f"{r['region']}|{r['collection_date']}\n")
            s = seqs[r["accession"]]
            for i in range(0, len(s), 60):
                handle.write(s[i:i + 60] + "\n")
    print(f"{lineage}: appended", flush=True)

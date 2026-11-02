"""Collect additional Alpha (B.1.1.7) isolates for the geographic task.

Geographic region and variant are badly confounded in GenBank: Omicron
records are overwhelmingly North American, Delta overwhelmingly European and
North American. Predicting region across all variants would therefore let a
model read the variant and infer the region from that, which is the circular
shortcut this project is built to avoid.

Holding the variant fixed removes the shortcut. Alpha has the widest regional
spread, so the geographic task uses Alpha alone. This script deepens the Alpha
sample so the task survives deduplication and quality filtering.

Run: python3 -I src/topup_alpha.py [cap_per_region] [max_pages]
"""
from __future__ import annotations

import collections
import csv
import sys
import time
import json
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import build_dataset as bd
import spike_reference as ref

CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 600
PAGES = int(sys.argv[2]) if len(sys.argv) > 2 else 90
LINEAGE = "B.1.1.7"

meta_path = Path("data/raw/metadata.csv")
fasta_path = Path("data/raw/spike_proteins.fasta")
existing = {r["accession"] for r in csv.DictReader(open(meta_path))}
print(f"{len(existing)} accessions already held", flush=True)

kept: list[dict] = []
per_region: collections.Counter = collections.Counter()
token, scanned = None, 0
for page in range(PAGES):
    params = {"filter.pangolin_classification": LINEAGE, "page_size": bd.PAGE_SIZE}
    if token:
        params["page_token"] = token
    payload = json.loads(bd.fetch(
        f"{bd.DATASETS}/virus/taxon/sars2/dataset_report?" + urllib.parse.urlencode(params)))
    reports = payload.get("reports") or []
    if not reports:
        break
    for rec in reports:
        scanned += 1
        iso, loc = rec.get("isolate") or {}, rec.get("location") or {}
        region, date = loc.get("geographic_region"), iso.get("collection_date")
        if not region or not date or rec.get("is_annotated") is False:
            continue
        if rec["accession"] in existing or per_region[region] >= CAP:
            continue
        kept.append({
            "accession": rec["accession"], "lineage": LINEAGE,
            "variant": ref.variant_from_pango(LINEAGE),
            "country": (loc.get("geographic_location") or "").split(":")[0].strip(),
            "region": region, "collection_date": date,
            "genome_length": rec.get("length"), "completeness": rec.get("completeness"),
        })
        per_region[region] += 1
        existing.add(rec["accession"])
    token = payload.get("next_page_token")
    if page % 10 == 0:
        print(f"  page {page:3d} scanned {scanned:6d} kept {len(kept):5d} {dict(per_region)}", flush=True)
    if not token:
        break
    time.sleep(bd.DATASETS_PAUSE)

print(f"\nscanned {scanned}, new records {len(kept)}: {dict(per_region)}", flush=True)
if not kept:
    sys.exit(0)

seqs = bd.fetch_spikes([r["accession"] for r in kept])
rows = [r for r in kept if r["accession"] in seqs]
print(f"{len(rows)} new records have a spike sequence", flush=True)

with open(meta_path, "a", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
    writer.writerows(rows)
with open(fasta_path, "a") as handle:
    for r in rows:
        handle.write(f">{r['accession']}|{r['lineage']}|{r['variant']}|"
                     f"{r['region']}|{r['collection_date']}\n")
        s = seqs[r["accession"]]
        for i in range(0, len(s), 60):
            handle.write(s[i:i + 60] + "\n")
print("appended", flush=True)

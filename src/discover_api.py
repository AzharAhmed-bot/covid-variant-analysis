"""Probe the NCBI Datasets v2 virus API and print its actual field names.

Run this first once outbound access to api.ncbi.nlm.nih.gov is available.
NCBI has renamed virus report fields between releases, so the mapping into the
canonical table is written only after seeing a real response.

Run: python3 -I src/discover_api.py
"""
import json
import sys

sys.path.insert(0, "src")
import acquire


def show(label, obj, limit=2500):
    print(f"\n{'=' * 70}\n{label}\n{'=' * 70}")
    print(json.dumps(obj, indent=2)[:limit])


try:
    raw = acquire._request(f"{acquire.DATASETS_BASE}/version")
    show("GET /version", json.loads(raw), 400)
except Exception as exc:
    print(f"version probe failed: {type(exc).__name__}: {exc}")

for lineage in ["B.1.1.7", "BA.1"]:
    try:
        records = acquire.api_dataset_report(lineage, page_size=2, max_records=2)
        print(f"\nlineage {lineage}: {len(records)} records")
        if records:
            show(f"first report record for {lineage}", records[0], 4000)
            print("\ntop-level keys:", sorted(records[0].keys()))
    except Exception as exc:
        print(f"{lineage} report failed: {type(exc).__name__}: {exc}")

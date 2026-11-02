"""Download the isolate-level dataset: spike proteins plus metadata."""
from __future__ import annotations

import collections
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import spike_reference as ref

UA = {"User-Agent": "strathmore-bioinformatics-coursework/1.0"}
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

LINEAGES = {
    "Alpha":   ["B.1.1.7"],
    "Delta":   ["B.1.617.2", "AY.4"],
    "Omicron": ["BA.1", "BA.2", "BA.5"],
}

PAGE_SIZE = 1000
MAX_PAGES_PER_LINEAGE = 25
CAP_PER_REGION = 150
EFETCH_BATCH = 100
EFETCH_PAUSE = 0.5
DATASETS_PAUSE = 0.25


def fetch(url: str, timeout: int = 300, retries: int = 4) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed after {retries} attempts: {url}") from last


def collect_metadata(lineage: str, target: int) -> list[dict]:
    """Page the dataset report, keeping usable records balanced across regions."""
    kept: list[dict] = []
    per_region: collections.Counter = collections.Counter()
    token, scanned = None, 0

    for _ in range(MAX_PAGES_PER_LINEAGE):
        params = {"filter.pangolin_classification": lineage, "page_size": PAGE_SIZE}
        if token:
            params["page_token"] = token
        payload = json.loads(
            fetch(f"{DATASETS}/virus/taxon/sars2/dataset_report?"
                  + urllib.parse.urlencode(params))
        )
        reports = payload.get("reports") or []
        if not reports:
            break

        for rec in reports:
            scanned += 1
            iso = rec.get("isolate") or {}
            loc = rec.get("location") or {}
            region = loc.get("geographic_region")
            date = iso.get("collection_date")
            if not region or not date or rec.get("is_annotated") is False:
                continue
            if per_region[region] >= CAP_PER_REGION:
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

        token = payload.get("next_page_token")
        if not token or len(kept) >= target:
            break
        time.sleep(DATASETS_PAUSE)

    print(f"  {lineage:12s} scanned {scanned:6d} -> kept {len(kept):5d}  "
          f"{dict(per_region)}", flush=True)
    return kept[:target]


def fetch_spikes(accessions: list[str]) -> dict[str, str]:
    """Return {accession: spike protein} using annotated CDS translations."""
    out: dict[str, str] = {}
    for i in range(0, len(accessions), EFETCH_BATCH):
        batch = accessions[i : i + EFETCH_BATCH]
        params = {"db": "nuccore", "id": ",".join(batch),
                  "rettype": "fasta_cds_aa", "retmode": "text"}
        text = fetch(f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(params),
                     timeout=600).decode("utf-8", "replace")
        cur, keep = None, False
        for line in text.splitlines():
            if line.startswith(">"):
                keep = "[protein=surface glycoprotein]" in line
                cur = line.split("|")[1].split("_prot")[0] if keep and "|" in line else None
                if cur:
                    out[cur] = []
            elif keep and cur:
                out[cur].append(line.strip())
        print(f"    efetch {i + len(batch):5d}/{len(accessions)}  "
              f"spikes so far {len(out)}", flush=True)
        time.sleep(EFETCH_PAUSE)
    return {k: "".join(v).replace("*", "").upper() for k, v in out.items()}


def main() -> None:
    target_per_variant = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    out_dir = Path("data/raw")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Collecting metadata", flush=True)
    records: list[dict] = []
    for variant, lineages in LINEAGES.items():
        per_lineage = max(1, target_per_variant // len(lineages))
        for lineage in lineages:
            records.extend(collect_metadata(lineage, per_lineage))

    print(f"\n{len(records)} records with metadata", flush=True)
    print("Fetching spike proteins", flush=True)
    seqs = fetch_spikes([r["accession"] for r in records])

    rows = [r for r in records if r["accession"] in seqs]
    print(f"\n{len(rows)} records have a spike sequence "
          f"({len(records) - len(rows)} without)", flush=True)

    meta_path = out_dir / "metadata.csv"
    with open(meta_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fasta_path = out_dir / "spike_proteins.fasta"
    with open(fasta_path, "w") as handle:
        for r in rows:
            handle.write(f">{r['accession']}|{r['lineage']}|{r['variant']}|"
                         f"{r['region']}|{r['collection_date']}\n")
            seq = seqs[r["accession"]]
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")

    counts = collections.Counter((r["variant"], r["region"]) for r in rows)
    print(f"\nwrote {meta_path} and {fasta_path}")
    print("\nvariant x region:")
    for (v, reg), n in sorted(counts.items()):
        print(f"  {v:8s} {reg:16s} {n}")


if __name__ == "__main__":
    main()

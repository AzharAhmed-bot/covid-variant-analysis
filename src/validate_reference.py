"""Sanity-check the reference spike translation against independent annotation."""
import sys
from collections import Counter

sys.path.insert(0, "src")
import spike_reference as ref

genome, clusters = sys.argv[1], sys.argv[2]
failures = []

spike = ref.load_reference_spike(genome)
print(f"Spike length      : {len(spike)} aa (expected {ref.EXPECTED_SPIKE_LENGTH})")
print(f"First 20 residues : {spike[:20]}")
print(f"Last 20 residues  : {spike[-20:]}")

known = {339: "G", 417: "K", 452: "L", 478: "T", 484: "E", 501: "N",
         570: "A", 614: "D", 681: "P", 716: "T", 982: "S", 1118: "D"}
print("\n[1] Known reference residues")
for pos, aa in sorted(known.items()):
    got = spike[pos - 1]
    if got != aa:
        failures.append(f"pos {pos}: expected {aa}, got {got}")
    print(f"  {'OK ' if got == aa else 'FAIL'} pos {pos:4d}: expected {aa}, got {got}")

all_clades = ref.load_covariants_mutations(clusters)

print("\n[2] Project clades: reference residues must match Wuhan-Hu-1")
for clade, variant in ref.CLADE_TO_VARIANT.items():
    muts = all_clades.get(clade, [])
    if not muts:
        failures.append(f"{clade}: no spike mutations found in CoVariants")
        print(f"  FAIL {clade:20s} not found")
        continue
    bad = [m for m in muts
           if m["pos"] <= len(spike) and len(m["ref"]) == 1
           and spike[m["pos"] - 1] != m["ref"]]
    n_rbd = sum(1 for m in muts if ref.in_rbd(m["pos"]))
    if bad:
        failures.append(f"{clade}: {bad}")
    print(f"  {'OK ' if not bad else 'FAIL'} {clade:20s} {variant:8s} "
          f"n_spike={len(muts):3d} n_RBD={n_rbd:3d} mismatches={len(bad)}")

print("\n[3] All clades: informational cross-check")
seen, agree, parent_annotated = set(), 0, []
for clade, muts in all_clades.items():
    for m in muts:
        if m["pos"] > len(spike) or len(m["ref"]) != 1 or not m["ref"].isalpha():
            continue
        key = (m["pos"], m["ref"])
        if key in seen:
            continue
        seen.add(key)
        if spike[m["pos"] - 1] == m["ref"]:
            agree += 1
        else:
            parent_annotated.append(f"{clade}: {m['ref']}{m['pos']}{m['alt']} "
                                    f"(Wuhan-Hu-1 has {spike[m['pos']-1]})")
total = agree + len(parent_annotated)
print(f"  {agree}/{total} unique (pos, ref) pairs match Wuhan-Hu-1 "
      f"({100*agree/total:.1f}%)")
print(f"  {len(parent_annotated)} annotated against a mutated parent lineage:")
for e in parent_annotated:
    print(f"    {e}")

covered = Counter()
for start, end in ref.SPIKE_DOMAINS.values():
    for p in range(start, end + 1):
        covered[p] += 1
gaps = [p for p in range(1, ref.EXPECTED_SPIKE_LENGTH + 1) if covered[p] == 0]
overlaps = [p for p, c in covered.items() if c > 1]
if gaps or overlaps:
    failures.append(f"domain partition: {len(gaps)} gaps, {len(overlaps)} overlaps")
print(f"\n[4] Domain partition over 1..{ref.EXPECTED_SPIKE_LENGTH}: "
      f"gaps={len(gaps)} overlaps={len(overlaps)}")
print(f"    RBD={ref.RBD} ({ref.RBD[1]-ref.RBD[0]+1} aa), "
      f"RBM={ref.RBM}, furin site={ref.FURIN_SITE}")

print("\n[5] Pango lineage mapping")
cases = {"B.1.1.7": "Alpha", "Q.1": "Alpha", "B.1.617.2": "Delta", "AY.4": "Delta",
         "AY.103": "Delta", "BA.1": "Omicron", "BA.5.2": "Omicron",
         "B.1.1.529": "Omicron", "B.1.177": None, "": None, None: None}
for lineage, expected in cases.items():
    got = ref.variant_from_pango(lineage)
    if got != expected:
        failures.append(f"variant_from_pango({lineage!r}) -> {got}, expected {expected}")
    print(f"  {'OK ' if got == expected else 'FAIL'} {str(lineage):12s} -> {got}")

print("\n" + "=" * 60)
if failures:
    print(f"CHECKS FAILED ({len(failures)})")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("ALL CHECKS PASSED")

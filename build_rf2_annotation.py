#!/usr/bin/env python3
"""
build_rf2_annotation.py

Convert the RF2-PPI (Zhang/Cong 2026 Science) .contacts files into the annotation
file format expected by annotation_test.py.

Input : best_models/<AC1>_<AC2>/<AC1>_S#__<AC2>_S#__<model>.contacts
        3 columns: residue_in_protein1, residue_in_protein2, AF2 contact probability
        (already filtered by the authors to distance <8 A and probability >0.6;
        residue numbers are relative to the FULL protein, not the segment)

Output: uniprot_id, aa_pos, annotation_id     with annotation_id = "<focal>_<partner>"
        BOTH directions are written (A->B and B->A), matching PIONEER's convention,
        so each protein's own interface residues are testable.

Multiple segment pairs for the same protein pair are merged.

Usage:
  python build_rf2_annotation.py --root .../contacts/best_models \
      --out rf2_annotation.tsv [--min-prob 0.6] [--rvas-mapped ASD_mapped.tsv]
"""
import argparse, os, glob, sys
from collections import defaultdict

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True, help="best_models directory")
    p.add_argument("--out", required=True)
    p.add_argument("--min-prob", type=float, default=0.6,
                   help="contact probability threshold (files are pre-filtered at 0.6)")
    p.add_argument("--rvas-mapped", default=None,
                   help="optional: restrict to proteins present in this variant table")
    a = p.parse_args()

    keep = None
    if a.rvas_mapped:
        import csv
        keep = {r["uniprot_id"] for r in csv.DictReader(open(a.rvas_mapped), delimiter="\t")}
        print(f"restricting to {len(keep):,} proteins with variants", file=sys.stderr)

    files = glob.glob(os.path.join(a.root, "*", "*.contacts"))
    print(f"{len(files):,} .contacts files", file=sys.stderr)

    # (focal, partner) -> set of residues
    iface = defaultdict(set)
    npair = skipped = nline = 0
    for i, f in enumerate(files):
        if i and i % 5000 == 0:
            print(f"  {i:,}/{len(files):,}", file=sys.stderr)
        base = os.path.basename(f)
        try:
            left, right = base.split("__")[0], base.split("__")[1]
            a1 = left.rsplit("_S", 1)[0]
            a2 = right.rsplit("_S", 1)[0]
        except Exception:
            skipped += 1
            continue
        if keep is not None and a1 not in keep and a2 not in keep:
            continue
        with open(f) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    r1, r2, pr = int(parts[0]), int(parts[1]), float(parts[2])
                except ValueError:
                    continue
                if pr < a.min_prob:
                    continue
                iface[(a1, a2)].add(r1)      # residues of protein 1 contacting protein 2
                iface[(a2, a1)].add(r2)      # and vice versa
                nline += 1
        npair += 1

    prots = {k[0] for k in iface}
    print(f"parsed {nline:,} contacts | {npair:,} segment-pair files | "
          f"{len(iface):,} directed interfaces | {len(prots):,} proteins", file=sys.stderr)
    if skipped:
        print(f"  skipped {skipped} unparseable filenames", file=sys.stderr)

    with open(a.out, "w") as o:
        o.write("uniprot_id\taa_pos\tannotation_id\n")
        n = 0
        for (focal, partner), residues in sorted(iface.items()):
            aid = f"{focal}_{partner}"
            for r in sorted(residues):
                o.write(f"{focal}\t{r}\t{aid}\n"); n += 1
    print(f"wrote {n:,} rows -> {a.out}", file=sys.stderr)

if __name__ == "__main__":
    main()

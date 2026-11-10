#!/usr/bin/env python3
"""Summarise Step 3 tblastn (VvECPG1 + chr01g05760 vs Vitaceae genomes).

For each genome: merge HSPs (either query, default run) into loci (same contig+strand, within 300 bp (VvECPG1 and 05760 are 811 bp apart)),
keep loci with best E <= 1e-3, extract the ORF that contains the best HSP (scan upstream to the nearest
in-frame Met with no stop between, downstream to the stop codon), and write
  tblastn/loci.tsv  (one row per locus; best E/bits per query, default and SEG; ORF length; Gly %)
  tblastn/orfs.fa   (translated ORFs, for re-scoring and alignment)
Single-exon ORF extraction is adequate here: both grape genes are single-exon (CDS within one exon).
"""
import csv, os, subprocess, sys, collections, shutil
ROOT = os.environ.get("VVECPG1_UPSTREAM_ROOT", ".")
D = f"{ROOT}/results/homology/tblastn"
GEN = f"{ROOT}/references/vitaceae_genomes"
BLAST = os.environ.get("BLAST_BIN", os.path.dirname(shutil.which("blastp") or "blastp")) or "."
SPECIES = ["Vitis_vinifera_PN40024", "Tetrastigma_loheri", "Ampelopsis_glandulosa",
           "Parthenocissus_tricuspidata", "Cissus_quadrangularis", "Arabidopsis_thaliana"]
QN = {"VvECPG1": "ECPG1", "Vitvi05_01chr01g05760": "c05760", "Vitvi05_01chr01g05760_t001": "c05760"}
CODE = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
RC = str.maketrans("ACGTacgtN", "TGCAtgcaN")

def read(path):
    if not os.path.exists(path): return []
    out = []
    for r in csv.reader(open(path), delimiter="\t"):
        out.append(dict(q=QN[r[0]], s=r[1], pid=float(r[2]), qs=int(r[6]), qe=int(r[7]),
                        ss=int(r[8]), se=int(r[9]), e=float(r[10]), bits=float(r[11]), frame=int(r[12])))
    return out

def fetch(sp, contig, a, b):
    seq = subprocess.run([f"{BLAST}/blastdbcmd", "-db", f"{GEN}/{sp}.fa", "-entry", contig,
                          "-range", f"{a}-{b}"], capture_output=True, text=True, check=True).stdout
    return "".join(seq.split("\n")[1:]).upper()

def slen(sp, contig):
    o = subprocess.run([f"{BLAST}/blastdbcmd", "-db", f"{GEN}/{sp}.fa", "-entry", contig, "-outfmt", "%l"],
                       capture_output=True, text=True, check=True).stdout
    return int(o.split()[0])

def orf(sp, h, flank=3000):
    lo, hi = min(h["ss"], h["se"]), max(h["ss"], h["se"])
    L = slen(sp, h["s"])
    a, b = max(1, lo - flank), min(L, hi + flank)
    s = fetch(sp, h["s"], a, b)
    if h["frame"] < 0:
        s = s.translate(RC)[::-1]; start = b - hi          # offset of HSP start in rc string
    else:
        start = lo - a
    end = start + (hi - lo + 1)
    # walk upstream codon by codon from HSP start to nearest Met before a stop
    i, met = start, None
    while i - 3 >= 0:
        i -= 3; c = CODE.get(s[i:i+3], "X")
        if c == "*": break
        if c == "M": met = i
    if CODE.get(s[start:start+3]) == "M" and met is None: met = start
    j = start
    while j + 3 <= len(s):
        c = CODE.get(s[j:j+3], "X")
        if c == "*": break
        j += 3
    st = met if met is not None else start
    prot = "".join(CODE.get(s[k:k+3], "X") for k in range(st, j, 3))
    complete = met is not None and j + 3 <= len(s)
    return prot, complete

rows, fa = [], []
for sp in SPECIES:
    hs = {t: read(f"{D}/{sp}.{t}.tsv") for t in ("default", "seg")}
    good = [h for h in hs["default"] if h["e"] <= 1e-3]
    good.sort(key=lambda h: (h["s"], h["frame"] > 0, min(h["ss"], h["se"])))
    loci = []
    for h in good:
        lo, hi = min(h["ss"], h["se"]), max(h["ss"], h["se"])
        if loci and loci[-1]["s"] == h["s"] and (loci[-1]["frame"] > 0) == (h["frame"] > 0) and lo - loci[-1]["hi"] < 300:
            loci[-1]["hi"] = max(loci[-1]["hi"], hi); loci[-1]["h"].append(h)
        else:
            loci.append(dict(s=h["s"], frame=h["frame"], lo=lo, hi=hi, h=[h]))
    for n, lc in enumerate(sorted(loci, key=lambda l: min(x["e"] for x in l["h"])), 1):
        rec = dict(species=sp, locus=f"{sp}_L{n}", contig=lc["s"], start=lc["lo"], end=lc["hi"],
                   strand="+" if lc["frame"] > 0 else "-")
        for t in ("default", "seg"):
            pool = lc["h"] if t == "default" else [x for x in hs["seg"] if x["s"] == lc["s"]
                     and (x["frame"] > 0) == (lc["frame"] > 0) and min(x["ss"], x["se"]) <= lc["hi"] + 100
                     and max(x["ss"], x["se"]) >= lc["lo"] - 100]
            for q in ("ECPG1", "c05760"):
                b = [x for x in pool if x["q"] == q]
                best = min(b, key=lambda x: x["e"]) if b else None
                rec[f"{q}_{t}_E"] = f"{best['e']:.1e}" if best else ""
                rec[f"{q}_{t}_bits"] = best["bits"] if best else 0
                rec[f"{q}_{t}_pid"] = best["pid"] if best else ""
                rec[f"{q}_{t}_qrange"] = f"{best['qs']}-{best['qe']}" if best else ""
        best = min(lc["h"], key=lambda x: x["e"])
        prot, complete = orf(sp, best)
        rec.update(orf_aa=len(prot), orf_complete=complete,
                   gly_pct=round(100 * prot.count("G") / max(1, len(prot)), 1),
                   closer_default="ECPG1" if rec["ECPG1_default_bits"] > rec["c05760_default_bits"] else "c05760",
                   closer_seg=("ECPG1" if rec["ECPG1_seg_bits"] > rec["c05760_seg_bits"] else "c05760")
                              if max(rec["ECPG1_seg_bits"], rec["c05760_seg_bits"]) > 0 else "none")
        rows.append(rec); fa.append(f">{rec['locus']} {lc['s']}:{lc['lo']}-{lc['hi']}({rec['strand']})\n{prot}\n")

with open(f"{D}/loci.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
open(f"{D}/orfs.fa", "w").writelines(fa)
c = collections.Counter(r["species"] for r in rows)
for sp in SPECIES:
    print(sp, c.get(sp, 0), "loci at E<=1e-3;",
          "SEG-surviving (E<=1e-3):",
          sum(1 for r in rows if r["species"] == sp and any(r[k] and float(r[k]) <= 1e-3 for k in ("ECPG1_seg_E", "c05760_seg_E"))))

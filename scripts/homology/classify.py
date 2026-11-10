#!/usr/bin/env python3
"""Classify VvECPG1 and relatives into the EA1-box protein classes of Uebler, Marton & Dresselhaus 2015
(Plant Reprod 28:183, doi:10.1007/s00497-015-0269-z):

  EAL  ~100 aa; classical N-terminal or internal signal peptide / N-terminal TM before the EA1-box; P-box; C-terminal
       Ala-rich A-box; grasses only.
  EAG  70-120 aa; EA1-box at the C-terminus with <= 8 aa after it, usually ending in K; glycine 35.6-59.5 % across the
       protein outside the box on both sides; no A-/P-box, no signal peptide, no TM.
  EAC  everything else (60->400 aa, no conserved feature but the box; many with signal peptide or TM).

Inputs (results/homology/ea1box_classes/): ea1box_set.fa (InterPro IPR039926 members of A. thaliana, maize, rice, grape
+ our T2T and Ampelopsis sequences), signalp6/prediction_results.txt, deeptmhmm/biolib_results/predicted_topologies.3line.
EA1-box coordinates = Pfam PF28138 envelope (hmmsearch --cut_ga). Paper labels are attached where Fig. 1 names the gene.
Outputs: features_and_classes.tsv, similarity_to_VvECPG1.tsv, README section printed to stdout.
"""
import csv, os, re, subprocess, sys, collections
from pathlib import Path
import argparse
a = argparse.ArgumentParser()
a.add_argument("--rerun-hmm", action="store_true", help="Requires HMMER 3.4 on PATH; predictors stay frozen")
a.add_argument("--features-only", action="store_true", help="Rebuild classification without Biopython pairwise alignment")
args = a.parse_args()
P = Path(__file__).resolve().parents[2]
ROOT = str(P)
O = str(P / "source_addenda/D_ea1box_classes")
OUT = P / "results/homology"
OUT.mkdir(parents=True, exist_ok=True)
HMM = os.environ.get("HMMER_BIN", "")
EA1 = f"{O}/EA1_PF28138.hmm"

# genes named in Uebler et al. 2015 Fig. 1 (gene ID -> paper name, class)
PAPER = {"Q5G8Z3": ("ZmEA1", "EAL"), "G1FF49": ("ZmEAL1", "EAL"), "B6TF57": ("ZmEAL2", "EAL"),
         "A0A0P0X8R5": ("OsEAL1", "EAL"), "Q6Z4F5": ("OsEAL2", "EAL"), "Q6Z4F7": ("OsEAL3", "EAL"),
         "Q9FLJ3": ("AtEAG1", "EAG"), "A0ACS0ZT27": ("AtEAG2a", "EAG"), "A0ACS0ZVV8": ("AtEAG2b", "EAG"),
         "A0A0N7KMP7": ("OsEAG1", "EAG"),
         "Q1G3K5": ("AtEAC1", "EAC"), "Q75IY9": ("OsEAC1", "EAC"), "Q0JN14": ("OsEAC2", "EAC"),
         "Q5ZB16": ("OsEAC3", "EAC"), "A0A0P0V2B0": ("OsEAC4", "EAC"), "Q2QN78": ("OsEAC5a", "EAC")}

seqs, order = {}, []
for l in open(f"{O}/ea1box_set.fa"):
    if l.startswith(">"):
        k = l[1:].strip(); order.append(k); seqs[k] = ""
    else:
        seqs[k] += l.strip()

# Pfam EA1 envelope per sequence
dom = f"{O}/ea1_hmmsearch.domtbl"
aln = f"{O}/ea1_hmmsearch.txt"
if args.rerun_hmm:
    dom = str(OUT / "ea1_hmmsearch.domtbl")
    aln = str(OUT / "ea1_hmmsearch.txt")
    executable = str(Path(HMM) / "hmmsearch") if HMM else "hmmsearch"
    subprocess.run([executable, "-E", "1", "--domE", "1", "--domtblout", dom, "-o", aln,
                    EA1, f"{O}/ea1box_set.fa"], check=True)
CORE = (48, 83)   # EA1-box proper within the 93-position PF28138 model (N-terminal hydrophobic part and A-box excluded)
box = {}
cur = None; blocks = collections.defaultdict(list)
for l in open(aln):
    m = re.match(r"^>> (\S+)", l)
    if m: cur = m.group(1); dn = -1; continue
    m = re.match(r"^\s+== domain (\d+)\s+score: ([\d.]+) bits;\s+conditional E-value: (\S+)", l)
    if m and cur: blocks[cur].append(dict(score=float(m.group(2)), ce=float(m.group(3)), h=[], s=[])); continue
    if cur and blocks[cur]:
        m = re.match(r"^\s+EA1\s+(\d+)\s+(\S+)\s+(\d+)\s*$", l)
        if m: blocks[cur][-1]["h"].append((int(m.group(1)), m.group(2))); continue
        m = re.match(r"^\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s*$", l)
        if m and cur.startswith(m.group(1)[:10]) and blocks[cur][-1]["h"] and len(blocks[cur][-1]["s"]) < len(blocks[cur][-1]["h"]):
            blocks[cur][-1]["s"].append((int(m.group(2)), m.group(3)))
for n, bl in blocks.items():
    b = max(bl, key=lambda x: x["score"])
    hp, sp_ = b["h"][0][0], b["s"][0][0]
    hseq = "".join(x[1] for x in b["h"]); sseq = "".join(x[1] for x in b["s"])
    hpos, spos, first, last = hp - 1, sp_ - 1, None, None
    for hc, sc in zip(hseq, sseq):
        if hc != ".": hpos += 1
        if sc != "-": spos += 1
        if hc != "." and CORE[0] <= hpos <= CORE[1] and sc != "-":
            first = spos if first is None else first; last = spos
    ga = any(float(x.split()[13]) >= 27 for x in open(dom) if not x.startswith("#") and x.split()[0] == n)
    box[n] = dict(score=b["score"], ievalue=b["ce"], ali=(first, last), above_GA=ga)
# SignalP 6
sp = {}
spf = f"{O}/signalp6/prediction_results.txt"
if os.path.exists(spf):
    for l in open(spf):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        sp[f[0].strip()] = dict(pred=f[1], p_other=float(f[2]), p_sp=float(f[3]), cs=f[-1] if len(f) > 4 else "")
# DeepTMHMM
tm = {}
lines = open(f"{O}/deeptmhmm/predicted_topologies.3line").read().split("\n")
for i in range(0, len(lines) - 2, 3):
    if lines[i].startswith(">"):
        name, typ = lines[i][1:].rsplit(" | ", 1)
        tm[name.strip()] = dict(type=typ.strip(), n_tm=len(re.findall(r"M+", lines[i + 2])),
                                sp_region=("S" in lines[i + 2]))

def gly(s): return 100 * s.count("G") / max(1, len(s))
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8,
          P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2, X=0)
def max_kd(s, w=19):
    """max mean Kyte-Doolittle hydropathy over a 19-aa window (TM-like if >= ~1.6); NaN if shorter than the window"""
    return round(max(sum(KD[c] for c in s[i:i + w]) / w for i in range(len(s) - w + 1)), 2) if len(s) >= w else float("nan")

rows = []
for k in order:
    s = seqs[k]; acc = k.split("|")[1]
    b = box.get(k)
    # EA1-box span: Pfam alignment start .. end (the HMM 48-83 core = the box proper)
    bs, be = (b["ali"] if b else (None, None))
    up, dn = (s[:bs - 1], s[be:]) if b else (s, "")
    tail = len(dn) if b else None
    spp = sp.get(k, {})
    has_sp = spp.get("pred", "").startswith("SP") or tm.get(k, {}).get("type") == "SP+TM" or tm.get(k, {}).get("type") == "SP"
    has_tm = tm.get(k, {}).get("n_tm", 0) > 0
    abox = bool(re.search(r"A{4,}|(?:A.?){5,}", dn)) if b else False
    # paper criteria
    if b is None:
        cls = "no EA1 box"
    elif (has_sp or has_tm) and abox and len(s) <= 140 and tail is not None:
        cls = "EAL"
    elif 70 <= len(s) <= 120 and tail <= 8 and not has_sp and not has_tm and gly(up) >= 35 and (gly(dn) >= 20 or tail <= 3):
        cls = "EAG"
    else:
        cls = "EAC"
    pname, pcls = PAPER.get(acc, ("", ""))
    rows.append(dict(id=k, paper_name=pname, paper_class=pcls, length=len(s), gly_pct_total=round(gly(s), 1),
                     ea1_box=f"{bs}-{be}" if b else "", ea1_ievalue=f"{b['ievalue']:.1e}" if b else "", pfam_above_GA=b["above_GA"] if b else False,
                     gly_pct_before_box=round(gly(up), 1), aa_after_box=tail, last_aa=s[-1],
                     ala_rich_tail=abox, signalp6=spp.get("pred", "NA"), signalp6_p_sp=spp.get("p_sp", ""),
                     deeptmhmm=tm.get(k, {}).get("type", "NA"), deeptmhmm_n_tm=tm.get(k, {}).get("n_tm", ""),
                     max_kd19_before_box=max_kd(up), class_by_criteria=cls))
with open(str(OUT / "features_and_classes.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

# check: do the criteria reproduce the paper's labels?
lab = [r for r in rows if r["paper_class"]]
agree = sum(r["paper_class"] == r["class_by_criteria"] for r in lab)
print(f"criteria vs paper labels: {agree}/{len(lab)} agree")
for r in lab:
    if r["paper_class"] != r["class_by_criteria"]:
        print("  disagree:", r["paper_name"], r["paper_class"], "->", r["class_by_criteria"])

# global similarity of VvECPG1 (and 05760) to every labelled representative (Biopython global, BLOSUM62, gap 10/0.5, free end gaps (= EMBOSS needle defaults))
if args.features_only:
    print("Feature reconstruction complete; pairwise similarity skipped (--features-only).")
    raise SystemExit(0)
from Bio import Align
from Bio.Align import substitution_matrices
AL = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=-10, extend_gap_score=-0.5, end_gap_score=0)  # EMBOSS needle defaults
B62 = substitution_matrices.load("BLOSUM62")
def needle(a, b):
    al = AL.align(seqs[a], seqs[b])[0]
    x, y = al[0], al[1]
    L = len(x); idn = sum(1 for p, q in zip(x, y) if p == q and p != "-")
    simc = sum(1 for p, q in zip(x, y) if p != "-" and q != "-" and B62[p][q] > 0)
    return idn, L, simc, al.score
qs = ["Vvi|T2T|VvECPG1_chr01g05750", "Vvi|T2T|chr01g05760"]
sim = []
for q in qs:
    for r in rows:
        if r["paper_class"] or r["id"].startswith("Vvi|") or r["id"].startswith("Agl|"):
            if r["id"] == q: continue
            i, L, s_, sc = needle(q, r["id"])
            # box-only identity: align the EA1-box spans
            sim.append(dict(query=q.split("|")[-1], target=r["paper_name"] or r["id"], target_class=r["paper_class"] or r["class_by_criteria"],
                            target_len=r["length"], global_identity_pct=round(100 * i / L, 1),
                            global_similarity_pct=round(100 * s_ / L, 1), needle_score=sc))
with open(str(OUT / "similarity_to_VvECPG1.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(sim[0]), delimiter="\t"); w.writeheader(); w.writerows(sim)

for q in ("VvECPG1_chr01g05750", "chr01g05760"):
    print(f"\n== {q}: best needle scores per class")
    for c in ("EAL", "EAG", "EAC"):
        cand = sorted([x for x in sim if x["query"] == q and x["target_class"] == c and not x["target"].startswith(("Vvi|", "Agl|"))],
                      key=lambda x: -x["needle_score"])[:3]
        print(" ", c, "; ".join(f"{x['target']} {x['needle_score']:.0f} ({x['global_identity_pct']}% id)" for x in cand))
print()
for r in rows:
    if r["id"].startswith(("Vvi|", "Agl|")) or r["paper_name"]:
        print("\t".join(str(r[k]) for k in ("id", "paper_name", "paper_class", "class_by_criteria", "length", "gly_pct_total",
                                            "gly_pct_before_box", "ea1_box", "aa_after_box", "last_aa", "ala_rich_tail",
                                            "signalp6", "signalp6_p_sp", "deeptmhmm", "deeptmhmm_n_tm", "max_kd19_before_box")))

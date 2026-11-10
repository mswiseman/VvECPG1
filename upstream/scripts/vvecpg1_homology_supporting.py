#!/usr/bin/env python3
"""VvECPG1 homology: the supporting steps that were first run interactively (2026-11-09/10), as one script.

  S1  CDS / protein of VvECPG1 (Vitvi05_01chr01g05750_t001) and Vitvi05_01chr01g05760_t001 from PN40024 T2T PN.fa +
      v5.1 GFF3 (single-exon CDS, + strand).
  S2  tblastn VvECPG1 vs the 12X.v2 assembly (legacy locus VIT_01s0011g04570) and vs the whole PN40024 T2T genome
      (finds the tandem paralog); blastp VvECPG1 vs chr01g05760.
  S3  Vitaceae (Vitis + Tetrastigma) hits from the nr BLASTP (Step 2): sequences pulled from the local nr and scored
      against VvECPG1 and chr01g05760, unmasked and SEG-masked (which grape gene do they resemble?).
  S4  Six-frame Pfam EA1 (PF28138) scan of +/-20 kb around the syntenic loci in Parthenocissus and Ampelopsis
      (is the ECPG1-type gene really absent in Parthenocissus?).
  S5  VvECPG1 / chr01g05760 vs maize ZmEA1 (UniProt Q5G8Z3): EA1 HMM coverage and blastp.

Needs: Step 2 output blastp_nr.default.tsv, Step 3 genome BLAST dbs (references/vitaceae_genomes, made by
vvecpg1_homology_step3_domains.sh) and Pfam-A.hmm (references/homology_dbs). Writes results/homology/supporting/.
Usage: envs/rnaseq/bin/python scripts/vvecpg1_homology_supporting.py
"""
import collections, csv, os, re, subprocess, shutil
ROOT = os.environ.get("VVECPG1_UPSTREAM_ROOT", ".")
R = f"{ROOT}/results/homology"
O = f"{R}/supporting"
B = os.environ.get("BLAST_BIN", os.path.dirname(shutil.which("blastp") or "blastp")) or "."
H = os.environ.get("HMMER_BIN", os.path.dirname(shutil.which("hmmsearch") or "hmmsearch")) or "."
GEN = f"{ROOT}/references/vitaceae_genomes"
NRDB = os.environ.get("BLASTDB", "references/homology_dbs")
LEGACY = os.environ.get("VVECPG1_LEGACY_FASTA", "references/12Xv2/grapevine_12Xv2.fa")
os.makedirs(O, exist_ok=True)
env = dict(os.environ, BLASTDB=NRDB)
CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, env=env, **kw).stdout
def translate(s): return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
def fasta(path):
    d, k = {}, None
    for l in open(path):
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = ""
        else: d[k] += l.strip()
    return d
log = open(f"{O}/summary.txt", "w")
def say(*a):
    print(*a); print(*a, file=log)

# ---------------- S1 ----------------
gff = f"{ROOT}/references/pn40024_v5.1/PN_T2T.v5.1.gff3"
cds = {}
for l in open(gff):
    f = l.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] == "CDS":
        m = re.search(r"Parent=([^;]+)", f[8])
        if m and m.group(1) in ("Vitvi05_01chr01g05750_t001", "Vitvi05_01chr01g05760_t001"):
            cds.setdefault(m.group(1), []).append((f[0], int(f[3]), int(f[4]), f[6]))
# PN.fa names chromosomes PN1..PN19 (= chr01..chr19 of the GFF3)
seqs = {}
for t, parts in cds.items():
    assert len(parts) == 1 and parts[0][3] == "+", parts   # both single-exon, + strand
    c, a, b, _ = parts[0]
    nt = "".join(run(["samtools", "faidx", f"{ROOT}/references/pn40024_v5.1/PN.fa", f"PN{int(c[3:])}:{a}-{b}"]).split("\n")[1:]).upper()
    p = translate(nt); assert p.endswith("*") and "*" not in p[:-1]
    seqs[t] = (c, a, b, nt, p[:-1])
e, f_ = seqs["Vitvi05_01chr01g05750_t001"], seqs["Vitvi05_01chr01g05760_t001"]
open(f"{O}/VvECPG1_cds.fa", "w").write(f">VvECPG1_CDS Vitvi05_01chr01g05750_t001 {e[0]}:{e[1]}-{e[2]}(+)\n{e[3]}\n")
open(f"{O}/VvECPG1_protein.fa", "w").write(f">VvECPG1 Vitvi05_01chr01g05750_t001 PN40024_T2T_v5.1 {e[0]}:{e[1]}-{e[2]}(+)\n{e[4]}\n")
open(f"{O}/Vitvi05_01chr01g05760_protein.fa", "w").write(f">Vitvi05_01chr01g05760_t001 neighbour {f_[0]}:{f_[1]}-{f_[2]}(+)\n{f_[4]}\n")
say(f"S1 VvECPG1 CDS {len(e[3])} bp = {len(e[3]) // 3} codons incl. stop -> protein {len(e[4])} aa, Gly {e[4].count('G')} "
    f"({100 * e[4].count('G') / len(e[4]):.1f} %); chr01g05760 {len(f_[4])} aa, Gly {100 * f_[4].count('G') / len(f_[4]):.1f} %")
for n in ("VvECPG1_cds.fa", "VvECPG1_protein.fa", "Vitvi05_01chr01g05760_protein.fa"):
    same = list(fasta(f"{O}/{n}").values()) == list(fasta(f"{R}/{n}").values())
    say(f"   {n} identical to the sequence used in the analysis: {same}")
Q = f"{R}/VvECPG1_protein.fa"; Q2 = f"{R}/Vitvi05_01chr01g05760_protein.fa"

# ---------------- S2 ----------------
F6 = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
ldb = f"{ROOT}/work/blastdb/grapevine_12Xv2"
if not os.path.exists(ldb + ".nsq") and not os.path.exists(ldb + ".00.nsq"):
    os.makedirs(os.path.dirname(ldb), exist_ok=True)
    run([f"{B}/makeblastdb", "-in", LEGACY, "-dbtype", "nucl", "-out", ldb])
leg = run([f"{B}/tblastn", "-query", Q, "-db", ldb, "-evalue", "1e-5", "-seg", "no", "-outfmt", F6, "-max_target_seqs", "20"])
open(f"{O}/tblastn_VvECPG1_vs_12Xv2.tsv", "w").write(leg)
top = leg.split("\n")[0].split("\t")
say(f"S2 legacy 12X.v2: top hit {top[1]}:{top[8]}-{top[9]} {top[2]} % identity over {top[3]} aa (q{top[6]}-{top[7]}), E {top[10]}")
t2t = run([f"{B}/tblastn", "-query", Q, "-db", f"{GEN}/Vitis_vinifera_PN40024.fa", "-evalue", "1e-3", "-seg", "no",
           "-outfmt", F6, "-max_target_seqs", "50"])
open(f"{O}/tblastn_VvECPG1_vs_PN40024_T2T.tsv", "w").write(t2t)
say(f"S2 whole T2T genome tblastn (E <= 1e-3): {len(t2t.strip().splitlines())} HSPs (expected: VvECPG1 itself + chr01g05760 CDS 4086170-4086448):")
for r in t2t.strip().split("\n"):
    x = r.split("\t"); say(f"   {x[1]}:{x[8]}-{x[9]} {x[2]} % over {x[3]} aa, E {x[10]}")
bp = run([f"{B}/blastp", "-query", Q, "-subject", Q2, "-seg", "no", "-outfmt", "0"])
open(f"{O}/blastp_VvECPG1_vs_chr01g05760.txt", "w").write(bp)
m = re.search(r"Expect = (\S+).*?Identities = (\d+)/(\d+).*?Positives = (\d+)/", bp, re.S)
say(f"S2 blastp VvECPG1 vs chr01g05760: {m.group(2)}/{m.group(3)} identities, {m.group(4)} positives, E {m.group(1)}")

# ---------------- S3 ----------------
ids = sorted({r[1].split("|")[1] for r in csv.reader(open(f"{R}/blastp_nr.default.tsv"), delimiter="\t")
              if r[15].split()[0] in ("Vitis", "Tetrastigma")})
open(f"{O}/vitaceae_ids.txt", "w").write("\n".join(ids) + "\n")
run([f"{B}/blastdbcmd", "-db", "nr", "-entry_batch", f"{O}/vitaceae_ids.txt", "-target_only", "-out", f"{O}/vitaceae_nr_hits.fa"])
run([f"{B}/makeblastdb", "-in", f"{O}/vitaceae_nr_hits.fa", "-dbtype", "prot", "-out", f"{O}/_vdb"])
qq = f"{O}/_q.fa"; open(qq, "w").write(open(Q).read() + open(Q2).read())
desc = {}
for l in open(f"{O}/vitaceae_nr_hits.fa"):
    if l.startswith(">"):
        a = l[1:].split()[0]; desc[a] = l[1:].strip().split(" ", 1)[1].split(" >")[0]
rows = []
for seg in ("no", "yes"):
    out = run([f"{B}/blastp", "-query", qq, "-db", f"{O}/_vdb", "-seg", seg, "-evalue", "1", "-max_target_seqs", "100",
               "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen"])
    b = collections.defaultdict(dict); sl = {}
    for l in out.strip().split("\n"):
        q, s, pid, L, qs, qe, ss, se, ev, bit, slen = l.split(); sl[s] = slen
        if s not in b[q] or float(bit) > b[q][s][0]: b[q][s] = (float(bit), pid, f"{qs}-{qe}")
    qn = [k for k in b if k != "VvECPG1"][0] if len(b) > 1 else "Vitvi05_01chr01g05760_t001"
    for s in sorted(sl):
        a = b["VvECPG1"].get(s, (0, "", "")); c = b.get(qn, {}).get(s, (0, "", ""))
        rows.append([s, desc.get(s, ""), sl[s], "seg" if seg == "yes" else "none", *a, *c,
                     "VvECPG1" if a[0] > c[0] else "chr01g05760"])
with open(f"{O}/vitaceae_ecpg1_vs_05760.tsv", "w") as f:
    f.write("subject\tdescription\tslen\tmasking\tECPG1_bits\tECPG1_pident\tECPG1_qrange\tc05760_bits\tc05760_pident\tc05760_qrange\tcloser_to\n")
    for r in rows: f.write("\t".join(map(str, r)) + "\n")
for p in os.listdir(O):
    if p.startswith(("_vdb", "_q")): os.remove(f"{O}/{p}")
tet = [r for r in rows if "Tetrastigma" in r[1] and r[3] == "none" and int(r[2]) <= 80]
say(f"S3 Vitaceae nr hits: {len(ids)} proteins; short Tetrastigma copies (<= 80 aa) closer to VvECPG1: "
    f"{sum(r[-1] == 'VvECPG1' for r in tet)}/{len(tet)} unmasked, "
    f"{sum(r[-1] == 'VvECPG1' for r in rows if 'Tetrastigma' in r[1] and r[3] == 'seg' and int(r[2]) <= 80)}/"
    f"{sum(1 for r in rows if 'Tetrastigma' in r[1] and r[3] == 'seg' and int(r[2]) <= 80)} SEG-masked")
same = open(f"{O}/vitaceae_ecpg1_vs_05760.tsv").read() == open(f"{R}/vitaceae_ecpg1_vs_05760.tsv").read()
say(f"   table identical to the one used in the analysis: {same}")

# ---------------- S4 ----------------
EA1 = f"{O}/EA1_PF28138.hmm"
open(EA1, "w").write(run([f"{H}/hmmfetch", f"{ROOT}/references/homology_dbs/Pfam-A.hmm", "PF28138.1"]))
with open(f"{O}/six_frame_EA1_scan.tsv", "w") as out:
    out.write("genome\tcontig\twindow\torf_id\tfull_E\tbits\n")
    for sp, contig, a, b in (("Parthenocissus_tricuspidata", "OZ252945.1", 4593000, 4634000),
                             ("Ampelopsis_glandulosa", "OZ253060.1", 4590000, 4633000)):
        s = "".join(run([f"{B}/blastdbcmd", "-db", f"{GEN}/{sp}.fa", "-entry", contig, "-range", f"{a}-{b}"]).split("\n")[1:]).upper()
        say(f"S4 {sp} {contig}:{a}-{b}: {len(s)} bp, N = {s.count('N')}")
        rc = s[::-1].translate(str.maketrans("ACGT", "TGCA"))
        orfs = []
        for st, seq in (("+", s), ("-", rc)):
            for fr in range(3):
                p = "".join(CODE.get(seq[i:i + 3], "X") for i in range(fr, len(seq) - 2, 3)); pos = 0
                for orf in p.split("*"):
                    if len(orf) >= 20:
                        nt = fr + 3 * pos; g = a + nt if st == "+" else a + len(s) - 1 - nt
                        orfs.append((f"{st}{fr}_{g}", orf))
                    pos += len(orf) + 1
        tf = f"{O}/_six.fa"; open(tf, "w").write("".join(f">{n}\n{o}\n" for n, o in orfs))
        tbl = f"{O}/_six.tbl"
        run([f"{H}/hmmsearch", "-E", "10", "--max", "--tblout", tbl, EA1, tf])
        hits = [l.split() for l in open(tbl) if not l.startswith("#")]
        for h in sorted(hits, key=lambda x: float(x[4])):
            out.write(f"{sp}\t{contig}\t{a}-{b}\t{h[0]}\t{h[4]}\t{h[5]}\n")
        strong = [h for h in hits if float(h[4]) < 1e-3]
        say(f"   six-frame ORFs >= 20 aa: {len(orfs)}; EA1 hits E < 1e-3: {len(strong)} "
            + "; ".join(f"{h[0]} E {float(h[4]):.1e}" for h in sorted(strong, key=lambda x: float(x[4]))))
        os.remove(tf); os.remove(tbl)

# ---------------- S5 ----------------
zm = "MSSCPAIVNMKDDDGIGAMGAAVAFAAMGVFGIYFLWPVVGPTSAGMMMKAPGAAGWVICRAVFEANPQLYFTILRTAGAAAAAATFAACSIAS"  # UniProt Q5G8Z3
zf = f"{O}/ZmEA1_Q5G8Z3.fa"; open(zf, "w").write(f">ZmEA1_Q5G8Z3\n{zm}\n")
q3 = f"{O}/_q3.fa"; open(q3, "w").write(open(Q).read() + open(Q2).read() + open(zf).read())
d = run([f"{H}/hmmsearch", "--cut_ga", "--domtblout", "/dev/stdout", "-o", "/dev/null", EA1, q3])
with open(f"{O}/EA1_coverage_vs_ZmEA1.tsv", "w") as f:
    f.write("protein\tbits\ti_evalue\thmm_from\thmm_to\tali_from\tali_to\n")
    for l in d.split("\n"):
        if l and not l.startswith("#"):
            x = l.split(); f.write(f"{x[0]}\t{x[13]}\t{x[12]}\t{x[15]}\t{x[16]}\t{x[17]}\t{x[18]}\n")
            say(f"S5 EA1 HMM: {x[0]} model {x[15]}-{x[16]} of 93, residues {x[17]}-{x[18]}, {x[13]} bits")
bz = run([f"{B}/blastp", "-query", Q, "-subject", zf, "-seg", "no", "-outfmt", "0"])
open(f"{O}/blastp_VvECPG1_vs_ZmEA1.txt", "w").write(bz)
m = re.search(r"Expect = (\S+).*?Identities = (\d+)/(\d+)", bz, re.S)
say(f"S5 blastp VvECPG1 vs ZmEA1: {m.group(2)}/{m.group(3)} identities, E {m.group(1)}")
os.remove(q3)
log.close()

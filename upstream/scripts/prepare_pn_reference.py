#!/usr/bin/env python3
"""Build the primary PN40024 T2T + v5.1 alignment genome and GTF from verified inputs.

Inputs (all checksum-verified downloads):
  PN.fa (Zenodo; PN1..PN19), PN_T2T_mit.fa/.gff3, PN_T2T_chl.fa, PN_T2T.v5.1.gff3,
  reports/reference/assembly_compare_*.pairs.tsv (exact-sequence identity tables).
Seqid maps are derived ONLY from exact sequence identity (PN<n> == T2T_ref chr<nn> ==
Ensembl CP1266xx). chr00 (unplaced, absent from PN.fa, redundant with chromosomes) is excluded.
Outputs in references/pn40024_v5.1/build/:
  genome_nuc_org.fa(.fai)  chr01..chr19 + MT + chl_<contig> (decoy, unannotated)
  genome_nuc.fa(.fai)      chr01..chr19 only (for organelle-effect test)
  pn_v51.gtf               exon/transcript records, gene_id = v5.1 gene ID; MT genes prefixed 'MT:'
  seqmap_gff_to_fasta.tsv, seqmap_ensembl_to_fasta.tsv, v51.validation.tsv, v51.genes.tsv
Usage: prepare_pn_reference.py --root <project_root>
"""
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd
import pysam


def write_fasta(records, path: Path) -> None:
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    tmp.replace(path)
    pysam.faidx(str(path))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    ref = root / "references/pn40024_v5.1"
    b = ref / "build"
    b.mkdir(exist_ok=True)
    rep = root / "reports/reference"

    # exact-identity maps
    g = pd.read_csv(rep / "assembly_compare_zenodo_vs_grapedia.pairs.tsv", sep="\t")
    g = g[(g.set_a == "zenodo_PN") & (g.set_b == "grapedia_T2T_ref") & (g.status == "identical_sequence")]
    pn_to_chr = dict(zip(g.seqid_a, g.seqid_b))
    assert len(pn_to_chr) == 19 and all(v.startswith("chr") for v in pn_to_chr.values()), pn_to_chr
    e = pd.read_csv(rep / "assembly_compare_zenodo_vs_ensembl.pairs.tsv", sep="\t")
    e = e[e.status == "identical_sequence"]
    assert len(e) == 19
    ens_map = pd.DataFrame({"ensembl_seqid": e.seqid_b, "fasta_seqid": e.seqid_a.map(pn_to_chr),
                            "evidence": "identical_normalized_sequence_md5"})
    ens_map.to_csv(b / "seqmap_ensembl_to_fasta.tsv", sep="\t", index=False)

    nuc = []
    with pysam.FastxFile(str(ref / "PN.fa")) as fh:
        for r in fh:
            nuc.append((pn_to_chr[r.name], r.sequence.upper()))
    nuc.sort()
    org = []
    with pysam.FastxFile(str(ref / "PN_T2T_mit.fa")) as fh:
        recs = list(fh)
        assert len(recs) == 1 and recs[0].name == "PN"
        org.append(("MT", recs[0].sequence.upper()))
    with pysam.FastxFile(str(ref / "PN_T2T_chl.fa")) as fh:
        for r in fh:
            org.append(("chl_" + r.name.replace(":", "_").replace("-", "_"), r.sequence.upper()))
    write_fasta(nuc, b / "genome_nuc.fa")
    write_fasta(nuc + org, b / "genome_nuc_org.fa")

    smap = pd.DataFrame({"gff_seqid": [c for _, c in sorted(pn_to_chr.items(), key=lambda x: x[1])],
                         "fasta_seqid": [c for _, c in sorted(pn_to_chr.items(), key=lambda x: x[1])]})
    smap.to_csv(b / "seqmap_gff_to_fasta.tsv", sep="\t", index=False)
    pd.DataFrame({"gff_seqid": ["PN"], "fasta_seqid": ["MT"]}).to_csv(b / "seqmap_mt.tsv", sep="\t", index=False)

    py = str(Path(__file__).resolve().parent / "gff_tools.py")
    run = lambda *args: subprocess.run(["python", py, *args], check=True)
    run("validate", "--gff", str(ref / "PN_T2T.v5.1.gff3"), "--fasta-fai", str(b / "genome_nuc_org.fa.fai"),
        "--seqmap", str(b / "seqmap_gff_to_fasta.tsv"), "--out-prefix", str(b / "v51"))
    run("validate", "--gff", str(ref / "PN_T2T_mit.gff3"), "--fasta-fai", str(b / "genome_nuc_org.fa.fai"),
        "--seqmap", str(b / "seqmap_mt.tsv"), "--out-prefix", str(b / "mt"))
    run("to-gtf", "--gff", str(ref / "PN_T2T.v5.1.gff3"), "--seqmap", str(b / "seqmap_gff_to_fasta.tsv"),
        "--keep-seqids-from", str(b / "genome_nuc_org.fa.fai"), "--source", "PN_T2T.v5.1",
        "--out", str(b / "v51_nuclear.gtf"))
    # MT GFF is GeSeq-style (gene records; CDS/exon without transcript parents): count by gene span.
    # Used only for organellar read accounting, never for DE interpretation.
    seen = {}
    with open(b / "mt.gtf.tmp", "w") as out:
        for line in open(ref / "PN_T2T_mit.gff3"):
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) != 9 or f[2] != "gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].strip(";").split(";") if "=" in kv)
            gid = "MT:" + attrs.get("ID", f"{f[3]}_{f[4]}")
            seen[gid] = seen.get(gid, 0) + 1
            if seen[gid] > 1:
                gid = f"{gid}.dup{seen[gid]}"
            a = f'gene_id "{gid}"; transcript_id "{gid}.t"; transcript_type "organellar_gene_span";'
            out.write(f"MT\tPN_T2T_mit\texon\t{f[3]}\t{f[4]}\t.\t{f[6]}\t.\t{a}\n")
    os.replace(b / "mt.gtf.tmp", b / "mt.gtf")
    with open(b / "pn_v51.gtf.tmp", "w") as out:
        for f in ("v51_nuclear.gtf", "mt.gtf"):
            out.write(open(b / f).read())
    os.replace(b / "pn_v51.gtf.tmp", b / "pn_v51.gtf")
    print("done")


if __name__ == "__main__":
    main()

#!/bin/bash
#SBATCH -p regular
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -J ecpg1s3
#SBATCH -o logs/homology/step3_domains.%j.log
# VvECPG1 homology, Step 3 + domain check (2026-11-10).
#  A. Domain search: current Pfam-A (hmmscan, gathering cut-offs AND permissive E<=10) and full NCBI CDD (rpsblast,
#     CD-Search default E<=0.01 AND permissive E<=10). Queries: VvECPG1, chr01g05760, two Tetrastigma copies.
#  B. tblastn of VvECPG1 + chr01g05760 vs Vitaceae genomes absent from nr as annotated proteins
#     (Ampelopsis, Parthenocissus, Cissus), with Tetrastigma loheri + PN40024 T2T as positive controls and
#     Arabidopsis TAIR10 as negative control. Default and SEG-masked.
# NB: no `| head` under pipefail (SIGPIPE killed job 1546).
set -euo pipefail
ROOT=${VVECPG1_UPSTREAM_ROOT:-.}
OUT=$ROOT/results/homology
REF=$ROOT/references/homology_dbs
GEN=$ROOT/references/vitaceae_genomes
BLAST=${BLAST_BIN:-$(dirname "$(command -v blastp)")}
HMM=${HMMER_BIN:-$(dirname "$(command -v hmmsearch)")}
DS=${DATASETS_BIN:-$(dirname "$(command -v datasets)")}
THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p $OUT/domains $OUT/tblastn $REF $GEN

# ---------- queries ----------
Q=$OUT/domains/queries.fa
cat $OUT/VvECPG1_protein.fa $OUT/Vitvi05_01chr01g05760_protein.fa > $Q
python3 - $OUT/vitaceae_nr_hits.fa >> $Q <<'EOF'
import sys
want={"KAN4971768.1","KAN5101018.1"}   # best T. harmandii / T. voinierianum full-length copies
keep=False
for l in open(sys.argv[1]):
    if l.startswith(">"):
        keep=l[1:].split()[0] in want
        if keep: l=">"+l[1:].split()[0]+"\n"
    if keep: sys.stdout.write(l)
EOF
grep -c ">" $Q

# ---------- A1. Pfam ----------
# Use the original Pfam 39.0 database, not current_release.
[ -s "$REF/Pfam-A.hmm.h3p" ] || { echo 'Restore and hmmpress the recorded Pfam 39.0 database first.' >&2; exit 2; }
zcat $REF/Pfam.version.gz > $OUT/domains/Pfam.version.txt
$HMM/hmmscan --cpu $THREADS --cut_ga --domtblout $OUT/domains/pfam_cutga.domtbl $REF/Pfam-A.hmm $Q > $OUT/domains/pfam_cutga.txt
$HMM/hmmscan --cpu $THREADS -E 10 --domE 10 --max --domtblout $OUT/domains/pfam_permissive.domtbl $REF/Pfam-A.hmm $Q > $OUT/domains/pfam_permissive.txt

# ---------- A2. CDD ----------
[ -s "$REF/cdd/Cdd.pal" ] || { echo 'Restore the recorded CDD database first.' >&2; exit 2; }
ls -l --time-style=long-iso $REF/cdd/cddid.tbl.gz > $OUT/domains/cdd.version.txt
F='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle'
$BLAST/rpsblast -query $Q -db $REF/cdd/Cdd -evalue 0.01 -seg no -comp_based_stats 1 -num_threads $THREADS -outfmt "$F" > $OUT/domains/cdd_e0.01.tsv
$BLAST/rpsblast -query $Q -db $REF/cdd/Cdd -evalue 10   -seg no -comp_based_stats 1 -num_threads $THREADS -outfmt "$F" > $OUT/domains/cdd_e10.tsv

# ---------- B. genomes ----------
declare -A ACC=( [Ampelopsis_glandulosa]=GCA_965260675.1 [Parthenocissus_tricuspidata]=GCA_965260695.1
                 [Cissus_quadrangularis]=GCA_002878655.1 [Tetrastigma_loheri]=GCA_056579005.1
                 [Arabidopsis_thaliana]=GCF_000001735.4 )
for sp in "${!ACC[@]}"; do
  a=${ACC[$sp]}
  if [ ! -s $GEN/$sp.fa ]; then
    $DS/datasets download genome accession $a --include genome --filename $GEN/$sp.zip --no-progressbar
    unzip -o -q $GEN/$sp.zip -d $GEN/tmp_$sp
    cat $GEN/tmp_$sp/ncbi_dataset/data/$a/*.fna > $GEN/$sp.fa
    rm -rf $GEN/tmp_$sp $GEN/$sp.zip
  fi
  [ -s $GEN/$sp.fa.ndb ] || $BLAST/makeblastdb -in $GEN/$sp.fa -dbtype nucl -parse_seqids -out $GEN/$sp.fa > /dev/null
done
ln -sf $ROOT/references/pn40024_v5.1/PN.fa $GEN/Vitis_vinifera_PN40024.fa
[ -s $GEN/Vitis_vinifera_PN40024.fa.ndb ] || $BLAST/makeblastdb -in $GEN/Vitis_vinifera_PN40024.fa -dbtype nucl -parse_seqids -out $GEN/Vitis_vinifera_PN40024.fa > /dev/null
for sp in "${!ACC[@]}"; do echo -e "$sp\t${ACC[$sp]}"; done > $OUT/tblastn/genomes.tsv
echo -e "Vitis_vinifera_PN40024\tPN40024 T2T v5.1 (local)" >> $OUT/tblastn/genomes.tsv

TQ=$OUT/tblastn/queries.fa
cat $OUT/VvECPG1_protein.fa $OUT/Vitvi05_01chr01g05760_protein.fa > $TQ
TF='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe slen'
for sp in Ampelopsis_glandulosa Parthenocissus_tricuspidata Cissus_quadrangularis Tetrastigma_loheri Arabidopsis_thaliana Vitis_vinifera_PN40024; do
  for seg in no yes; do
    tag=$([ $seg = yes ] && echo seg || echo default)
    $BLAST/tblastn -query $TQ -db $GEN/$sp.fa -seg $seg -evalue 10 -max_target_seqs 500 -num_threads $THREADS \
      -outfmt "$TF" > $OUT/tblastn/${sp}.${tag}.tsv
  done
  echo "$sp done: $(wc -l < $OUT/tblastn/${sp}.default.tsv) / $(wc -l < $OUT/tblastn/${sp}.seg.tsv) HSPs"
done
echo ALL DONE

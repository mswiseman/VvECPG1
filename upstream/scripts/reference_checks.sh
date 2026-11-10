#!/usr/bin/env bash
# Reference reconciliation + build (idempotent). Reproduces reports/reference/* and references/*/build.
# Usage: scripts/reference_checks.sh <project_root>
set -euo pipefail
ROOT=$1; cd "$ROOT"
E=$ROOT/envs/rnaseq/bin
PY=$E/python
R=references; OUT=reports/reference; mkdir -p "$OUT" references/star "$R/ensembl_r62/build"
log() { echo "[$(date -u +%FT%TZ)] $*" >&2; }

# 1. Grapedia archives -> read-only extracted copies
mkdir -p $R/grapedia/extracted
for z in $R/grapedia/*.zip; do unzip -n -q "$z" -d $R/grapedia/extracted; done
chmod 0440 $R/grapedia/extracted/*

# 2. exact sequence identity: Zenodo PN.fa vs Ensembl r62 vs Grapedia T2T_ref (+ organelles)
[[ -s $OUT/assembly_compare_zenodo_vs_ensembl.pairs.tsv ]] || $PY scripts/compare_assemblies.py \
  --set zenodo_PN=$R/pn40024_v5.1/PN.fa \
  --set ensembl_r62=$R/ensembl_r62/Vitis_vinifera.ASM3070453v1.dna.toplevel.fa.gz \
  --out-prefix $OUT/assembly_compare_zenodo_vs_ensembl
[[ -s $OUT/assembly_compare_zenodo_vs_grapedia.pairs.tsv ]] || $PY scripts/compare_assemblies.py \
  --set zenodo_PN=$R/pn40024_v5.1/PN.fa --set grapedia_T2T_ref=$R/grapedia/extracted/T2T_ref.fasta \
  --set zenodo_chl=$R/pn40024_v5.1/PN_T2T_chl.fa --set zenodo_mit=$R/pn40024_v5.1/PN_T2T_mit.fa \
  --out-prefix $OUT/assembly_compare_zenodo_vs_grapedia

# 3. chr00 (present only in T2T_ref) redundancy vs chromosomes
if [[ ! -s $OUT/chr00_vs_PN.paf ]]; then
  $E/samtools faidx $R/grapedia/extracted/T2T_ref.fasta chr00 > $OUT/chr00.fa.tmp
  $E/minimap2 -c -t 16 -x asm5 --secondary=no $R/pn40024_v5.1/PN.fa $OUT/chr00.fa.tmp > $OUT/chr00_vs_PN.paf.tmp 2>/dev/null
  mv $OUT/chr00_vs_PN.paf.tmp $OUT/chr00_vs_PN.paf; rm -f $OUT/chr00.fa.tmp
fi
$PY - "$OUT/chr00_vs_PN.paf" > $OUT/chr00_redundancy.tsv <<'EOF'
import sys, collections
L = 9748653
cov = bytearray(L); bins = collections.Counter()
for l in open(sys.argv[1]):
    f = l.rstrip("\n").split("\t")
    de = [float(x[5:]) for x in f[12:] if x.startswith("de:f:")][0]   # requires minimap2 -c
    qs, qe = int(f[2]), int(f[3])
    cov[qs:qe] = b"\x01" * (qe - qs)
    b = ">=99.9" if de <= 0.001 else "99-99.9" if de <= 0.01 else "98-99" if de <= 0.02 else "95-98" if de <= 0.05 else "<95"
    bins[b] += qe - qs
print("metric\tbp\tfraction_of_chr00")
print(f"chr00_length\t{L}\t1")
print(f"covered_by_primary_alignments_to_chr01_19\t{sum(cov)}\t{sum(cov)/L:.4f}")
for k in [">=99.9", "99-99.9", "98-99", "95-98", "<95"]:
    print(f"aligned_bp_gap_compressed_identity_{k}\t{bins[k]}\t{bins[k]/L:.4f}")
EOF

# 4. primary genome + GTF + GFF validation + CDS integrity
B=$R/pn40024_v5.1/build
[[ -s $B/pn_v51.gtf ]] || $PY scripts/prepare_pn_reference.py --root "$ROOT"
[[ -s $B/v51.cds_integrity.tsv ]] || $PY scripts/check_cds_integrity.py --gff $R/pn40024_v5.1/PN_T2T.v5.1.gff3 \
  --fasta $B/genome_nuc_org.fa --seqmap $B/seqmap_gff_to_fasta.tsv --out $B/v51.cds_integrity.tsv
[[ -s $B/pn_v51_ensembl_seqids.gtf ]] || awk -F'\t' 'NR==FNR{if(FNR>1)m[$2]=$1; next} ($1 in m){$1=m[$1]; print}' \
  OFS='\t' $B/seqmap_ensembl_to_fasta.tsv $B/v51_nuclear.gtf > $B/pn_v51_ensembl_seqids.gtf
[[ -s $R/ensembl_r62/build/ensembl_r62.gtf ]] || zcat $R/ensembl_r62/Vitis_vinifera.ASM3070453v1.62.gtf.gz > $R/ensembl_r62/build/ensembl_r62.gtf

# 4b. chr00-inclusive sensitivity genome (+ Grapedia T2T_ref chr00 and its 882 v5.1 genes)
if [[ ! -s $B/pn_v51_chr00.gtf ]]; then
  ( cat $B/genome_nuc_org.fa; $E/samtools faidx $R/grapedia/extracted/T2T_ref.fasta chr00 ) > $B/genome_nuc_org_chr00.fa.tmp
  mv $B/genome_nuc_org_chr00.fa.tmp $B/genome_nuc_org_chr00.fa; $E/samtools faidx $B/genome_nuc_org_chr00.fa
  (cat $B/seqmap_gff_to_fasta.tsv; printf 'chr00\tchr00\n') > $B/seqmap_gff_to_fasta_chr00.tsv
  $PY scripts/gff_tools.py to-gtf --gff $R/pn40024_v5.1/PN_T2T.v5.1.gff3 --seqmap $B/seqmap_gff_to_fasta_chr00.tsv \
    --keep-seqids-from $B/genome_nuc_org_chr00.fa.fai --source PN_T2T.v5.1 --out $B/v51_with_chr00.gtf
  cat $B/v51_with_chr00.gtf $B/mt.gtf > $B/pn_v51_chr00.gtf
fi

# 5. STAR indexes (through SLURM when available)
run() { if command -v srun >/dev/null; then srun -p regular -c 16 --mem=64G "$@"; else "$@"; fi; }
run scripts/build_star_index.sh $B/genome_nuc_org.fa $B/pn_v51.gtf references/star/pn_v51 16 93 &
run scripts/build_star_index.sh $B/genome_nuc.fa $B/v51_nuclear.gtf references/star/pn_v51_nuc 16 93 &
run scripts/build_star_index.sh $R/ensembl_r62/Vitis_vinifera.ASM3070453v1.dna.toplevel.fa.gz \
  $R/ensembl_r62/Vitis_vinifera.ASM3070453v1.62.gtf.gz references/star/ensembl_r62 16 93 &
run scripts/build_star_index.sh $B/genome_nuc_org_chr00.fa $B/pn_v51_chr00.gtf references/star/pn_v51_chr00 16 93 &
wait
log "reference checks complete"

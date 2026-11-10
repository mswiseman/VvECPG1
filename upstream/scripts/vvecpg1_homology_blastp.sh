#!/usr/bin/env bash
# VvECPG1 homology search, step 2: BLASTP against annotated protein databases.
# Usage: sbatch scripts/vvecpg1_homology_blastp.sh
#
# Query: VvECPG1 (75 aa, Gly 17.3%), translated from PN40024 T2T v5.1 CDS chr01:4,084,956-4,085,183(+).
# The protein is short and compositionally biased, so every search is run twice:
#   .default  composition-based statistics on (BLAST+ default, -comp_based_stats 2), no low-complexity masking
#   .seg      the same, plus SEG masking of the query (-seg yes)
# A hit that disappears under SEG masking was driven by the biased segment, not by shared sequence.
# E-value is deliberately loose (10) so marginal hits are visible and can be judged, not silently dropped.
#SBATCH -p regular
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -J ecpg1bp
#SBATCH -o logs/homology/blastp.%j.log
set -euo pipefail
# NB: never pipe into `head` here - SIGPIPE + pipefail kills the job (exit 13); use sed -n.

ROOT=${VVECPG1_UPSTREAM_ROOT:-.}
OUT=$ROOT/results/homology
BLAST=${BLAST_BIN:-$(dirname "$(command -v blastp)")}
: "${BLASTDB:?Set BLASTDB to the archived 2026-04-07 nr/refseq database directory}"
export BLASTDB
VIRIDIPLANTAE=33090   # NCBI taxonomy id; -taxids descends the subtree in BLAST v5 databases
THREADS=${SLURM_CPUS_PER_TASK:-8}
FMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp staxids sscinames sskingdoms stitle'

mkdir -p "$OUT" "$ROOT/logs/homology"
QUERY=$OUT/VvECPG1_protein.fa
[ -s "$QUERY" ] || { echo "missing query $QUERY (write it first)" >&2; exit 1; }

echo "=== database versions ==="
for db in nr refseq_protein; do
  echo "--- $db"; $BLAST/blastdbcmd -db $db -info 2>&1 | sed -n 1,6p || true
done

run () {  # run <db> <tag> <extra blastp args...>
  local db=$1 tag=$2; shift 2
  echo "=== blastp $db ($tag) ==="
  $BLAST/blastp -query "$QUERY" -db "$db" -taxids $VIRIDIPLANTAE \
      -evalue 10 -max_target_seqs 1000 -num_threads "$THREADS" \
      -outfmt "$FMT" "$@" > "$OUT/blastp_${db}.${tag}.tsv"
  # pairwise alignments for the top hits, so alignment extent can be inspected by eye
  $BLAST/blastp -query "$QUERY" -db "$db" -taxids $VIRIDIPLANTAE \
      -evalue 10 -max_target_seqs 50 -num_threads "$THREADS" \
      -outfmt 0 "$@" > "$OUT/blastp_${db}.${tag}.aln.txt"
  echo "$db $tag: $(wc -l < "$OUT/blastp_${db}.${tag}.tsv") hit lines"
}

for db in refseq_protein nr; do
  run "$db" default
  run "$db" seg -seg yes
done

echo "=== done $(date -Is) ==="

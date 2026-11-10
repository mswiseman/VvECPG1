#!/usr/bin/env bash
# Quantify one SRA run with Salmon for the focal-gene atlas (SLURM array task).
# Usage: sbatch --array=1-N%6 scripts/sra_atlas_quant.sh <run_list.txt>
# Per run: first 10 M spots via fastq-dump, read 1 only (same treatment for SE and PE studies),
# Salmon selective alignment vs PN40024 v5.1 transcripts with whole-genome decoys, libType A.
# Keeps quant.sf.gz + aux meta_info.json; FASTQ deleted. Skips runs already done.
set -euo pipefail
ROOT="${VVECPG1_UPSTREAM_ROOT:-$PWD}"
LIST=$1
RUN=$(sed -n "${SLURM_ARRAY_TASK_ID:-${RUN_INDEX:-1}}p" "$LIST")
OUT=$ROOT/results/sra_atlas/quant/$RUN
[ -s "$OUT/quant.sf.gz" ] && { echo "$RUN done already"; exit 0; }
TMP="$ROOT/work/sra_atlas/$RUN"
mkdir -p "$TMP" "$OUT"
cd "$TMP"
for try in 1 2 3; do
  if "${FASTQ_DUMP:-fastq-dump}" -X 10000000 --split-files --skip-technical -O "$TMP" "$RUN"; then break; fi
  echo "fastq-dump attempt $try failed"; rm -f "$TMP"/*.fastq; sleep 60
done
R1=$(ls "$TMP"/${RUN}_1.fastq 2>/dev/null || ls "$TMP"/${RUN}.fastq)
"${SALMON:-salmon}" quant -i $ROOT/references/salmon_v51/index_k31 -l A -r "$R1" \
  -p 8 --validateMappings --gcBias -o "$TMP/salmon" --quiet
gzip -c "$TMP/salmon/quant.sf" > "$OUT/quant.sf.gz"
cp "$TMP/salmon/aux_info/meta_info.json" "$OUT/meta_info.json"
cp "$TMP/salmon/lib_format_counts.json" "$OUT/lib_format_counts.json"
rm -rf "$TMP"
echo "$RUN ok"

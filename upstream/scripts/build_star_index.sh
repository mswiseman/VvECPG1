#!/usr/bin/env bash
# Build a STAR genome index atomically (writes to <out>.tmp, then renames).
# Usage: build_star_index.sh <genome.fa[.gz]> <annotation.gtf[.gz]|none> <out_dir> <threads> <sjdbOverhang>
set -euo pipefail
fa=$1; gtf=$2; out=$3; threads=$4; overhang=${5:-93}
ENV=$(cd "$(dirname "$0")/.." && pwd)/envs/rnaseq/bin
if [[ -s $out/SA && -s $out/genomeParameters.txt ]]; then echo "index exists: $out"; exit 0; fi
tmp=$out.tmp; rm -rf "$tmp"; mkdir -p "$tmp"
work=$(mktemp -d "${out}.inputs.XXXXXX"); trap 'rm -rf "$work"' EXIT
case $fa in *.gz) zcat "$fa" > "$work/genome.fa"; fa=$work/genome.fa;; esac
args=()
if [[ $gtf != none ]]; then
  case $gtf in *.gz) zcat "$gtf" > "$work/ann.gtf"; gtf=$work/ann.gtf;; esac
  args=(--sjdbGTFfile "$gtf" --sjdbOverhang "$overhang")
fi
"$ENV/STAR" --version
"$ENV/STAR" --runMode genomeGenerate --runThreadN "$threads" --genomeDir "$tmp" \
  --genomeFastaFiles "$fa" --genomeSAindexNbases 14 --outTmpDir "$work/_STARtmp" \
  --outFileNamePrefix "$tmp/" "${args[@]}"
rm -rf "$out"; mv "$tmp" "$out"
echo "built $out"

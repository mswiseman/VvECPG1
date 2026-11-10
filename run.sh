#!/usr/bin/env bash
# Run from any directory. Output is confined to this package.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/.cache
export XDG_CACHE_HOME="$PWD/results/.cache"
export MPLCONFIGDIR="$PWD/results/.cache/matplotlib"
mkdir -p "$MPLCONFIGDIR"
PYTHON="${PYTHON:-python3}"
RSCRIPT="${RSCRIPT:-Rscript}"
case "${1:-check}" in
  check|check-data|restore-data) ;;
  *)
    if [[ ! -d inputs || ! -d reference_results || ! -d source_addenda ]]; then
      echo 'Restore the full release data first: bash run.sh restore-data /path/to/VvECPG1_submission_code_complete_2026-11-10.zip' >&2
      exit 2
    fi
    ;;
esac
case "${1:-check}" in
  check)
    "$PYTHON" scripts/check_repository.py
    "$RSCRIPT" scripts/check_environment.R
    ;;
  check-data)
    "$PYTHON" scripts/check_inputs.py
    ;;
  restore-data)
    "$PYTHON" scripts/restore_data.py "${2:?Supply the path to the full release ZIP}"
    ;;
  verify-reference)
    "$PYTHON" scripts/check_inputs.py
    export VVECPG1_RESULTS=reference_results
    "$RSCRIPT" scripts/coexpression/verify.R
    "$RSCRIPT" scripts/functional/verify.R
    "$PYTHON" scripts/hormone/summarize.py
    ;;
  figures-reference)
    export VVECPG1_RESULTS=reference_results
    "$RSCRIPT" scripts/figures/atlas.R
    "$RSCRIPT" scripts/figures/combined_S4.R
    ;;
  reproduce)
    "$PYTHON" scripts/check_inputs.py
    "$RSCRIPT" scripts/check_environment.R
    export VVECPG1_RESULTS=results
    "$RSCRIPT" scripts/coexpression/run.R
    "$RSCRIPT" scripts/coexpression/verify.R
    "$RSCRIPT" scripts/functional/run.R
    "$RSCRIPT" scripts/functional/verify.R
    "$PYTHON" scripts/hormone/summarize.py
    "$RSCRIPT" scripts/figures/atlas.R
    "$RSCRIPT" scripts/figures/combined_S4.R
    ;;
  paralog-reference)
    export VVECPG1_RESULTS=reference_results
    "$RSCRIPT" scripts/paralog/correlations.R
    "$PYTHON" scripts/paralog/expression.py
    ;;
  hormone-figure)
    "$PYTHON" scripts/hormone/summarize.py
    "$PYTHON" scripts/hormone/plot.py
    ;;
  atlas-rebuild)
    "$PYTHON" scripts/atlas/rebuild.py --matrix
    ;;
  homology-features)
    "$PYTHON" scripts/homology/classify.py --features-only
    ;;
  homology-similarity)
    "$PYTHON" scripts/homology/classify.py
    ;;
  prepare-upstream)
    "$PYTHON" scripts/prepare_upstream.py
    ;;
  programme-rebuild)
    "$PYTHON" scripts/prepare_upstream.py
    VVECPG1_UPSTREAM_ROOT="$PWD/results/upstream_workspace" "$RSCRIPT" upstream/scripts/focal_24h_programme.R
    ;;
  *) echo 'Usage: bash run.sh {check|check-data|restore-data|verify-reference|figures-reference|reproduce|paralog-reference|hormone-figure|atlas-rebuild|homology-features|homology-similarity|prepare-upstream|programme-rebuild}' >&2; exit 2 ;;
esac

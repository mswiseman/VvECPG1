#!/usr/bin/env bash
# Extract zip members atomically (CRC-checked by unzip), skipping already-extracted members.
# Originals are never modified. Extracted files are made read-only.
# Usage: scripts/extract_inputs.sh <zip> <dest_dir> [parallel_jobs] [member_glob]
set -euo pipefail
zip=$1; dest=$2; jobs=${3:-6}; glob=${4:-*}
mkdir -p "$dest/.tmp"
extract_one() {
  local zip=$1 dest=$2 m=$3 base
  base=$(basename "$m")
  [[ -s $dest/$base ]] && return 0
  local tmp; tmp=$(mktemp -d "$dest/.tmp/x.XXXXXX")
  if unzip -q -j "$zip" "$m" -d "$tmp"; then
    mv "$tmp/$base" "$dest/$base"; chmod a-w "$dest/$base"; rmdir "$tmp"
    echo "ok $base"
  else
    echo "FAILED $m" >&2; rm -rf "$tmp"; return 1
  fi
}
export -f extract_one
unzip -Z1 "$zip" | grep -v '/$' | grep -E "$(printf '%s' "$glob" | sed 's/\./\\./g; s/\*/.*/g')" \
  | xargs -P "$jobs" -I{} bash -c 'extract_one "$@"' _ "$zip" "$dest" {}
rmdir "$dest/.tmp" 2>/dev/null || true
echo "done: $(ls "$dest" | wc -l) files in $dest"

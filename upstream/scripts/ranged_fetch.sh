#!/usr/bin/env bash
# Parallel byte-range download (for throttled hosts such as Zenodo). Atomic; caller verifies checksum.
# Usage: ranged_fetch.sh <url> <dest> <size_bytes> [n_chunks]
set -euo pipefail
url=$1; dest=$2; size=$3; n=${4:-16}
tmp=$(mktemp -d "${dest}.chunks.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
chunk=$(( (size + n - 1) / n ))
pids=()
for ((i=0; i<n; i++)); do
  s=$(( i * chunk )); e=$(( s + chunk - 1 )); (( e >= size )) && e=$(( size - 1 ))
  (( s > e )) && break
  curl -fsSL --retry 8 --retry-delay 5 -r "$s-$e" -o "$tmp/$(printf %04d $i)" "$url" & pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
cat "$tmp"/* > "$dest.part"
[[ $(stat -c %s "$dest.part") -eq $size ]] || { echo "size mismatch for $dest" >&2; exit 1; }
mv "$dest.part" "$dest"

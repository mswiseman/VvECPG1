#!/usr/bin/env bash
# Download and verify reference resources (idempotent, atomic, checksum-verified).
# Usage: scripts/download_references.sh [project_root]
# Writes references/<set>/<file> plus references/<set>/SOURCES.tsv
# (url, file, size, md5, sha256, verified_against, retrieved_utc).
set -euo pipefail
ROOT=${1:-$(cd "$(dirname "$0")/.." && pwd)}
cd "$ROOT"

ZENODO_RECORD=13741133
ZENODO_DOI="10.5281/zenodo.${ZENODO_RECORD}"
ZENODO_FILES=(PN.fa PN_T2T.v5.1.gff3 PN40024.cds.fa PN_protein.fa correelation.list.txt
              PN_T2T_chl.fa PN_T2T_chl.gff3 PN_T2T_mit.fa PN_T2T_mit.gff3)
ENS=https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-62
ENS_FILES=(
  "fasta/vitis_vinifera/dna|Vitis_vinifera.ASM3070453v1.dna.toplevel.fa.gz"
  "gtf/vitis_vinifera|Vitis_vinifera.ASM3070453v1.62.gtf.gz"
  "gff3/vitis_vinifera|Vitis_vinifera.ASM3070453v1.62.gff3.gz"
  "fasta/vitis_vinifera/cdna|Vitis_vinifera.ASM3070453v1.cdna.all.fa.gz"
  "fasta/vitis_vinifera/pep|Vitis_vinifera.ASM3070453v1.pep.all.fa.gz"
)
GRAPEDIA=https://grapedia.org/wp-content/uploads
GRAPEDIA_FILES=(2024/11/all_on_T2T_ref_equivalences_filtered_no_copies_symbols.csv
                2024/11/all_on_T2T_ref_equivalences_filtered_symbols.csv
                2023/11/T2T_ref.zip
                2024/11/PN40024_5.1_on_T2T_ref_with_names.zip
                2024/07/5.1_on_T2T_ref_gene_coordinates.zip)

log() { echo "[$(date -u +%FT%TZ)] $*" >&2; }

# record <dir> <url> <file> <verified_against>
record() {
  local dir=$1 url=$2 f=$3 ver=$4 tsv=$1/SOURCES.tsv
  [[ -s $tsv ]] || printf 'url\tfile\tsize_bytes\tmd5\tsha256\tverified_against\tretrieved_utc\n' > "$tsv"
  grep -qP "\t${f}\t" "$tsv" && return 0
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$url" "$f" "$(stat -c %s "$dir/$f")" \
    "$(md5sum "$dir/$f" | cut -d' ' -f1)" "$(sha256sum "$dir/$f" | cut -d' ' -f1)" \
    "$ver" "$(date -u +%FT%TZ)" >> "$tsv"
}

# fetch <url> <dest>: atomic download (resume-safe; skips if dest exists)
fetch() {
  local url=$1 dest=$2
  if [[ -s $dest ]]; then log "exists: $dest"; return 0; fi
  log "GET $url"
  curl -fsSL --retry 5 --retry-delay 10 -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}

# --- Zenodo (PN40024 T2T + v5.1) ---
d=references/pn40024_v5.1; mkdir -p "$d"
meta=$d/zenodo_record_${ZENODO_RECORD}.json
[[ -s $meta ]] || { curl -fsSL "https://zenodo.org/api/records/${ZENODO_RECORD}" > "$meta.part" && mv "$meta.part" "$meta"; }
for f in "${ZENODO_FILES[@]}"; do
  url="https://zenodo.org/records/${ZENODO_RECORD}/files/${f}?download=1"
  size=$(python3 -c "import json;d=json.load(open('$meta'));print([x['size'] for x in d['files'] if x['key']=='$f'][0])")
  if [[ ! -s $d/$f ]]; then log "GET (ranged) $url"; "$(dirname "$0")/ranged_fetch.sh" "$url" "$d/$f" "$size" 16; fi
  want=$(python3 -c "import json,sys;d=json.load(open('$meta'));print([x['checksum'] for x in d['files'] if x['key']=='$f'][0].split(':')[1])")
  got=$(md5sum "$d/$f" | cut -d' ' -f1)
  [[ $want == "$got" ]] || { log "MD5 MISMATCH $f want=$want got=$got"; mv "$d/$f" "$d/$f.bad"; exit 1; }
  record "$d" "$url" "$f" "zenodo_api_md5 (DOI ${ZENODO_DOI})"
  chmod a-w "$d/$f"
done

# --- Ensembl Plants release 62 (vendor reference lineage, GCA_030704535.1) ---
d=references/ensembl_r62; mkdir -p "$d"
for spec in "${ENS_FILES[@]}"; do
  sub=${spec%%|*}; f=${spec##*|}
  fetch "$ENS/$sub/$f" "$d/$f"
  cks=$d/CHECKSUMS.$(echo "$sub" | tr '/' '_')
  [[ -s $cks ]] || curl -fsSL "$ENS/$sub/CHECKSUMS" > "$cks"
  want=$(awk -v f="$f" '$3==f{print $1+0" "$2+0}' "$cks")
  got=$(sum "$d/$f" | awk '{print $1+0" "$2+0}')
  [[ -n $want && $want == "$got" ]] || { log "SUM MISMATCH $f want='$want' got='$got'"; mv "$d/$f" "$d/$f.bad"; exit 1; }
  record "$d" "$ENS/$sub/$f" "$f" "ensembl_CHECKSUMS_bsd_sum"
  chmod a-w "$d/$f"
done

# --- Grapedia ID equivalences (no published checksum; record our hashes) ---
d=references/grapedia; mkdir -p "$d"
for rel in "${GRAPEDIA_FILES[@]}"; do
  f=$(basename "$rel")
  fetch "$GRAPEDIA/$rel" "$d/$f"
  record "$d" "$GRAPEDIA/$rel" "$f" "none_published(recorded_sha256)"
  chmod a-w "$d/$f"
done
log "all references downloaded and verified"

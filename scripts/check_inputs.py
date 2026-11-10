#!/usr/bin/env python3
"""Validate bundled data/checkpoint hashes without modifying scientific inputs."""
from pathlib import Path
import csv
import hashlib

ROOT = Path(__file__).resolve().parents[1]
with (ROOT / 'provenance/input_manifest.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
errors = []
for row in rows:
    path = ROOT / row['packaged_file']
    if not path.is_file():
        errors.append(f'Missing: {row["packaged_file"]}')
    elif hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
        errors.append(f'Changed: {row["packaged_file"]}')
if errors:
    raise SystemExit('Analysis inputs are missing or changed. Restore the complete release ZIP first:\n  python3 scripts/restore_data.py /path/to/VvECPG1_submission_code_complete_2026-11-10.zip\n\n' + '\n'.join(errors[:20]) + (f'\n... {len(errors)} errors total' if len(errors)>20 else ''))
print(f'PASS: {len(rows)} input/reference files match their SHA-256 hashes.')

#!/usr/bin/env python3
"""Check tracked repository files independently of optional analysis data."""
from pathlib import Path
import hashlib
P=Path(__file__).resolve().parents[1]
errors=[];n=0
for line in (P/'SHA256SUMS').read_text().splitlines():
    expected,name=line.split(maxsplit=1);f=P/name;n+=1
    if not f.is_file():errors.append('Missing repository file: '+name)
    elif hashlib.file_digest(f.open('rb'),'sha256').hexdigest()!=expected:errors.append('Changed repository file: '+name)
if errors:raise SystemExit('\n'.join(errors))
print(f'PASS: {n} repository files match SHA-256 checksums. Analysis data are checked separately.')

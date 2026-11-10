#!/usr/bin/env python3
"""Restore only analysis inputs from the checksum-pinned full release ZIP.

Never overwrites repository code or existing files with different contents.
Python 3.11+; standard library only. No network access or shell execution.
"""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json,os,shutil,stat,tempfile,zipfile
P=Path(__file__).resolve().parents[1]
ALLOWED={'inputs','reference_results','source_addenda'}
def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def allowed_path(name):
    path=PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or '\\' in name:
        raise ValueError('Unsafe archive path: '+name)
    parts=path.parts
    if len(parts)<3 or parts[1] not in ALLOWED:return None
    return Path(*parts[1:])
def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('archive',type=Path);args=a.parse_args()
    archive=args.archive.expanduser().resolve();record=json.loads((P/'data_archive.json').read_text())
    if archive.stat().st_size!=record['bytes']:raise SystemExit('Archive size differs from the pinned full package.')
    print('Checking full-archive SHA-256...',flush=True)
    if digest(archive)!=record['sha256']:raise SystemExit('Archive SHA-256 does not match data_archive.json; nothing extracted.')
    with (P/'provenance/input_manifest.tsv').open() as f:
        import csv
        expected={r['packaged_file']:r['sha256'] for r in csv.DictReader(f,delimiter='\t')}
    with zipfile.ZipFile(archive) as z:
        members={}
        for info in z.infolist():
            rel=allowed_path(info.filename)
            if rel is None or info.is_dir():continue
            if stat.S_ISLNK(info.external_attr>>16):raise SystemExit('Archive symlink refused: '+info.filename)
            key=rel.as_posix()
            if key in members:raise SystemExit('Duplicate archive path: '+key)
            members[key]=info
        if set(members)!=set(expected):raise SystemExit('Archive input roster differs from the pinned input manifest.')
        # Check all existing destinations before writing anything. Re-running is safe.
        for name in members:
            dst=P/name
            if any(x.is_symlink() for x in [dst,*dst.parents] if x==P or P in x.parents):
                raise SystemExit('Destination symlink refused: '+name)
            if dst.exists() and (not dst.is_file() or digest(dst)!=expected[name]):
                raise SystemExit('Refusing to replace an existing different file: '+name)
        written=0
        for name,info in members.items():
            dst=P/name
            if dst.exists():continue
            dst.parent.mkdir(parents=True,exist_ok=True)
            fd,tmp=tempfile.mkstemp(prefix='restore-',suffix='.partial',dir=dst.parent)
            tmp=Path(tmp)
            try:
                with os.fdopen(fd,'wb') as out,z.open(info) as src:shutil.copyfileobj(src,out)
                if digest(tmp)!=expected[name]:raise SystemExit('Input SHA-256 mismatch: '+name)
                tmp.replace(dst);written+=1
            finally:
                tmp.unlink(missing_ok=True)
    print(f'PASS: {len(expected)} input/source files verified; {written} restored. Repository code unchanged.')
if __name__=='__main__':main()

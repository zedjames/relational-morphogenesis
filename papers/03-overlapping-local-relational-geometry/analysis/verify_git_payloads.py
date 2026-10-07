"""Publication-only Git-blob gate; numerical verification itself does not require Git."""
import argparse
import hashlib
import io
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', action='store_true', help='Validate staged objects before committing')
    args = parser.parse_args()
    repo = ROOT.parents[1]
    records = json.loads((ROOT/'verification/source_export_manifest.json').read_text())['files']
    prefix = ROOT.relative_to(repo).as_posix()+'/'
    revision = '' if args.index else 'HEAD'
    requests = ''.join(revision+':'+prefix+rel+'\n' for rel in records).encode()
    data = subprocess.check_output(['git', '-C', str(repo), 'cat-file', '--batch'], input=requests)
    stream = io.BytesIO(data)
    for rel, record in records.items():
        header = stream.readline().decode().strip().split()
        assert len(header) == 3 and header[1] == 'blob', (rel, header)
        content = stream.read(int(header[2]))
        assert stream.read(1) == b'\n'
        assert not content.startswith(b'version https://git-lfs.github.com/spec/v1'), 'Committed LFS pointer: '+rel
        assert len(content) == record['bytes'] and hashlib.sha256(content).hexdigest() == record['sha256'], rel
    assert not stream.read()
    print(f'PASS {len(records)} actual Git payloads match frozen byte counts/SHA-256; zero committed LFS pointers')

if __name__ == '__main__':
    main()

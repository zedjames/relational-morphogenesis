"""Frozen source mapping, complete release hashes and publication boundary."""
import gzip
import hashlib
import json
import re
import zipfile
from pathlib import Path
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

def verify():
    source_path = ROOT / 'verification/frozen_final08_artifact_manifest.json'
    source = json.loads(source_path.read_text())['files']
    export = json.loads((ROOT / 'verification/source_export_manifest.json').read_text())
    assert sha(source_path) == export['source_manifest_sha256']
    assert export['source_files'] == len(source) == 473
    assert export['source_bytes'] == sum(r['bytes'] for r in source.values())
    seen = set()
    for rel, rec in export['files'].items():
        original = rec['source_path'].split('/reviewer_release/', 1)[1]
        assert original not in seen and source[original] == {k: rec[k] for k in ('bytes', 'sha256')}
        seen.add(original)
        path = ROOT / rel
        assert not Path(rel).is_absolute() and '..' not in Path(rel).parts
        assert not path.is_symlink() and path.resolve().is_relative_to(ROOT)
        assert path.is_file() and path.stat().st_size == rec['bytes'] and sha(path) == rec['sha256'], rel
        data = path.read_bytes()
        assert not data.startswith(b'version https://git-lfs.github.com/spec/v1'), rel
        if path.suffix == '.npz':
            with zipfile.ZipFile(path) as archive:
                assert archive.testzip() is None and archive.namelist()
                assert all(n.endswith('.npy') for n in archive.namelist())
        elif path.suffix == '.gz':
            gzip.decompress(data)
    for rec in export['excluded_files']:
        rel = rec['source_path']
        assert rel not in seen and source[rel] == {k: rec[k] for k in ('bytes', 'sha256')}
        assert rel.startswith(('manuscript/', 'rmmo_paper3_final08/publication/manuscript/')) or rel == 'rmmo_paper3_referee05/.gitattributes'
        seen.add(rel)
    assert seen == set(source)
    assert export['exported_files'] == len(export['files'])
    assert export['exported_bytes'] == sum(r['bytes'] for r in export['files'].values())
    manifest_path = ROOT / 'verification/release_manifest.json'
    release = json.loads(manifest_path.read_text())['files']
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
              and '__pycache__' not in p.parts
              and not p.relative_to(ROOT).as_posix().startswith('reproducibility/rmmo_paper3_final08/completion/')}
    assert actual == set(release) | {'verification/release_manifest.json'}, 'Unmanifested or missing release files'
    for rel, rec in release.items():
        path = ROOT / rel
        assert not Path(rel).is_absolute() and '..' not in Path(rel).parts
        assert not path.is_symlink() and path.resolve().is_relative_to(ROOT)
        assert path.stat().st_size == rec['bytes'] and sha(path) == rec['sha256'], rel
    findings = []
    rules = [r'/Us' + r'ers/[^\s"\']+', r'/pri' + r'vate/tmp/[^\s"\']+',
             r'github\.com[/:]zedjames/qpci-planetary',
             r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',
             r'\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{40,}|AKIA[A-Z0-9]{16})\b']
    for path in ROOT.rglob('*'):
        rel = path.relative_to(ROOT)
        assert not path.is_symlink(), str(rel)
        assert not any(p in ('.git', '.aws', '.ssh', '.env') for p in rel.parts), str(rel)
        assert path.suffix.lower() not in ('.h5td', '.h5', '.hdf5'), str(rel)
        if path.is_file():
            data = path.read_bytes()
            texts = [data.decode('utf-8', errors='ignore')]
            if path.suffix == '.gz':
                texts.append(gzip.decompress(data).decode('utf-8', errors='ignore'))
            elif path.suffix == '.npz':
                with zipfile.ZipFile(path) as archive:
                    texts.extend(archive.read(n).decode('utf-8', errors='ignore') for n in archive.namelist())
            elif path.suffix == '.pdf':
                reader = PdfReader(path)
                texts.append(str(reader.metadata))
                texts.extend(page.extract_text() or '' for page in reader.pages)
            for rule in rules:
                if any(re.search(rule, text) for text in texts):
                    findings.append(str(rel))
    assert not findings, findings
    print(f"PASS frozen source mapping: {len(export['files'])} files / {export['exported_bytes']} bytes; boundary and checksums")
    return export

if __name__ == '__main__':
    verify()

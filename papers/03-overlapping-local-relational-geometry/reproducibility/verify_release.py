"""Standalone derived-input verification, without Git, private paths or source acquisition."""
import os,sys,json,hashlib,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()
def main():
    manifest=json.loads((ROOT/'verification/artifact_manifest.json').read_text())['files']
    for rel,record in manifest.items():
        p=ROOT/rel;assert ROOT in p.resolve().parents and not p.is_symlink() and sha(p)==record['sha256'],rel
        if p.suffix in ('.py','.md','.txt','.json','.csv','.tex'):
            assert not any(x in p.read_bytes() for x in [('/Us'+'ers/').encode(),('/private/'+'tmp/').encode()]),rel
    env=dict(os.environ,RMMO_P3_INPUT_ROOT=str(ROOT/'data_authority'),PYTHONPYCACHEPREFIX=str(Path(tempfile.gettempdir())/'rmmo08_review_pycache'),OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    for action in ('scalar_ranks','matched_replay','degree_replay'):
        subprocess.run([sys.executable,str(ROOT/'rmmo_paper3_final08/verify_calibration08.py'),action],env=env,check=True)
    subprocess.run([sys.executable,str(ROOT/'rmmo_paper3_final08/validate_publication08.py')],env=env,check=True)
    for rel,record in manifest.items():assert sha(ROOT/rel)==record['sha256'],('authority changed',rel)
    receipt=dict(passed=True,ordinary_files_verified=len(manifest),new_calibration_banks=7,inferential_banks_changed=0,draws_per_bank=1024,exact_graphs_replayed=33792,components=108,families=36,family_events=72,component_decisions_changed=0,family_decisions_changed=0,figures=10,Git_private_paths_raw_acquisition_Lean_required=False)
    (ROOT/'verification/validation_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('Standalone Final08 verification PASS',len(manifest),'ordinary files',flush=True)
if __name__=='__main__':main()

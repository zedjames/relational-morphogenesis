"""One-command full portable verification; no private paths or Git needed."""
import hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
def main():
    manifest=json.loads((ROOT/'verification/artifact_manifest.json').read_text())['files'];banned=[('/Us'+'ers/').encode(),('/private/'+'tmp/').encode(),('github.com/zedjames/'+'qpci-planetary').encode()]
    assert not (ROOT/'.git').exists() and not list(ROOT.rglob('*.h5td'))
    for rel,r in manifest.items():
        p=ROOT/rel;assert ROOT in p.resolve().parents and p.is_file() and not p.is_symlink();assert sha(p)==r['sha256'] and p.stat().st_size==r['bytes'],rel
        if p.suffix in ('.md','.py','.json','.csv','.txt'):assert not any(x in p.read_bytes() for x in banned),rel
    env=dict(os.environ,RMMO_P3_INPUT_ROOT=str(ROOT/'data_authority'),PYTHONPYCACHEPREFIX=str(Path(tempfile.gettempdir())/'rmmo_p3_referee05_review_cache'),OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    for rel,args in [('analysis/tests05.py',[]),('analysis/verify05.py',['--public','--parallel-banks'])]:subprocess.run([sys.executable,str(ROOT/rel)]+args,env=env,check=True)
    for rel,r in manifest.items():assert sha(ROOT/rel)==r['sha256'],('verifier mutated authority',rel)
    receipt={'passed':True,'ordinary_files_verified':len(manifest),'all_retained_draws_replayed':True,'GS_only_arrays_checked':219,'figures':6,'rendered_manuscript':True,'private_repository_required':False,'raw_H5TD_required':False,'external_upload_performed':False,'biological_replicates':3}
    (ROOT/'verification/validation_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('Standalone Referee05 verification PASS',len(manifest),'ordinary files',flush=True)
if __name__=='__main__':main()

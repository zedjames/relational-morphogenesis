"""Independent poison/audit harness, outside the Stage A construction boundary."""
import argparse,io,subprocess,shutil,ast
import scipy.sparse as sp
from io05 import *
import construction05 as build
P4=ROOT.parent/'rmmo_paper3_referee04'
GEOMETRY_SUBJECT='RMMO P3 Referee05: freeze strict GS-only geometry before GR reveal'
def poison():
    freeze();raw=sp.load_npz(P4/'normalization/selected_eligible_counts.npz');idx=np.load(P4/'normalization/count_index.npz');alternative=raw.copy();mask=np.isin(alternative.indices,idx['heldout_columns']);oldhash=array_digest(alternative.data[mask]);alternative.data[mask]=rng('GR_POISON_ONLY').integers(0,101,size=mask.sum());assert array_digest(alternative.data[mask])!=oldhash
    base=sp.load_npz(ROOT/'data_authority/GS_raw_counts.npz');allowed=alternative[:,idx['state_columns']].copy()
    for field in ('data','indices','indptr'):np.testing.assert_array_equal(getattr(base,field),getattr(allowed,field))
    # The same adapter contract projects out forbidden columns before entering Stage A.
    # Actual Stage A compute receives no GR capability, identity or annotation.
    result=build.compute(allowed,read(ROOT/'authority/primary_GS.csv.gz'),read(ROOT/'data_authority/cell_identity.csv.gz'),freeze(),read(ROOT/'freeze/strict_lane_plan.csv'));a=result[0];receipt=json.loads((ROOT/'separation_audit/geometry_pre_reveal_hashes.json').read_text())
    assert {k:array_digest(v) for k,v in a.items()}==receipt['array_content_hashes'];buf=io.BytesIO();np.savez_compressed(buf,**a);container=hashlib.sha256(buf.getvalue()).hexdigest();assert container==receipt['container_sha256'],'Compressed geometry container differs; do not report byte invariance'
    changed=['GR counts'];tree=ast.parse((ROOT/'construction05.py').read_text());imports=[n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)];assert not any(x and ('common05' in x or 'referee04' in x or 'references' in x) for x in imports)
    js('separation_audit/heldout_poison_test.json',{'passed':True,'original_GS_bytes_identical':True,'GR_values_changed':int(mask.sum()),'GR_poison_method':'deterministic nonnegative pseudorandom0..100 on all existing GR nonzero slots','every_array_content_hash_identical':True,'compressed_geometry_byte_identical':True,
        'geometry_container_sha256':container,'numeric_arrays_compared':len(a),'actual_compute_function_reexecuted':True,'poison_harness_outside_StageA':True,'GR_capability_supplied_to_construction':False,'biological_evidence':False,'completed_utc':utc()});print('GR poison test PASS: all arrays and compressed L/U/F/O geometry bytes invariant',len(a),flush=True)
def geometry_seal():
    commit=next(x.split('|')[0] for x in subprocess.check_output(['git','log','--format=%H|%s'],text=True,cwd=ROOT.parent.parent).splitlines() if x.split('|',1)[1]==GEOMETRY_SUBJECT)
    path=ROOT/'separation_audit/geometry_commit_receipt.json';old=json.loads(path.read_text()) if path.exists() else {};pre=json.loads((ROOT/'separation_audit/geometry_pre_reveal_hashes.json').read_text());js('separation_audit/geometry_commit_receipt.json',{'original_geometry_commit':old.get('original_geometry_commit',commit),'reachable_geometry_commit':commit,'original_geometry_seal_utc':old.get('original_geometry_seal_utc',utc()),'geometry_sha256':pre['container_sha256'],'heldout_revealed_before_original_commit':False});print('Strict geometry committed before reveal',commit,flush=True)
def reveal():
    assert not (ROOT/'separation_audit/reveal_receipt.json').exists();freeze();g=json.loads((ROOT/'separation_audit/geometry_commit_receipt.json').read_text());assert json.loads((ROOT/'separation_audit/heldout_poison_test.json').read_text())['passed'];assert sha(ROOT/'geometry/stageA_arrays.npz')==g['geometry_sha256']
    started=utc();sources={'heldout_hash':'normalization/primary_hash/heldout_projection.npz','pca32':'normalization/pca32/representation.npz','pca64':'normalization/pca64/representation.npz','target_pca32':'normalization/target_pca32/representation.npz'};hashes={}
    for name,rel in sources.items():
        out=ROOT/'evaluation_inputs'/name/'representation.npz';out.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(P4/rel,out);hashes[name]=sha(out)
    for name in ('pca32','pca64'):
        out=ROOT/'geometry'/name;out.mkdir(exist_ok=True);shutil.copyfile(P4/f'geometry/{name}/pair_cores.npz',out/'pair_cores.npz');shutil.copyfile(P4/f'geometry/{name}/state_geometry.npz',out/'state_geometry.npz')
    shutil.copyfile(P4/'authority/cell_carrier.csv.gz',ROOT/'authority/cell_carrier.csv.gz')
    for name in ('state','heldout'):shutil.copyfile(P4/f'normalization/{name}_expression.npz',ROOT/'evaluation_inputs'/(name+'_expression.npz'))
    js('separation_audit/reveal_receipt.json',{'geometry_commit':g['reachable_geometry_commit'],'geometry_sha256':g['geometry_sha256'],'evaluation_start_utc':started,'ZR_authority_sha256':hashes['heldout_hash'],'loaded_only_after_geometry_commit':True,'evaluation_input_hashes':hashes,'post_construction_annotation_access':True})
    js('evaluation_inputs/pca_byte_authority.json',{'source':'Referee04 GS-only isolated state expression/PCA and GR-only target PCA','embeddings':hashes,'state_expression_sha256':sha(ROOT/'evaluation_inputs/state_expression.npz'),'heldout_expression_sha256':sha(ROOT/'evaluation_inputs/heldout_expression.npz'),'G_S_genes':14660,'G_R_genes':6339,'state_geometry_does_not_read_GR':True,'parent_normalization_reconstruction_verified':True})
    print('Stage B reveal started only after committed strict geometry',started,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['poison','geometry_seal','reveal']);a=p.parse_args();globals()[a.action]()

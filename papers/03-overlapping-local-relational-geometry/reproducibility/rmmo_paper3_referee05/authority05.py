"""Pre-construction input extraction and independent identity audit, not Stage A."""
import argparse,shutil,subprocess
import scipy.sparse as sp
from io05 import *
P4=ROOT.parent/'rmmo_paper3_referee04'
def seal():
    subject='RMMO P3 Referee05: freeze GS-only strict construction before outcomes';commit=next(x.split('|')[0] for x in subprocess.check_output(['git','log','--format=%H|%s'],text=True,cwd=ROOT.parent.parent).splitlines() if x.split('|',1)[1]==subject)
    p=ROOT/'freeze/freeze_receipt.json';old=json.loads(p.read_text()) if p.exists() else {}
    js('freeze/freeze_receipt.json',{'original_plan_commit':old.get('original_plan_commit',commit),'reachable_plan_commit':commit,'original_seal_utc':old.get('original_seal_utc',utc()),
       'plan_sha256':{p:sha(ROOT/p) for p in ('freeze/analysis_plan.json','freeze/analysis_plan.md','freeze/strict_lane_plan.csv')},'wording':'The strict feature-separation construction and evaluation plan was committed before execution of the newly specified reconstruction.'});print('Strict plan sealed',commit,flush=True)
def prepare(directive,manuscript):
    freeze();assert not (ROOT/'authority/parent_authority.json').exists()
    parents=[f'relational_morphogenesis{i}' for i in range(12,18)]+[f'rmmo_paper3_referee0{i}' for i in range(1,5)]
    files={str(p.relative_to(ROOT.parent)):{'sha256':sha(p),'bytes':p.stat().st_size} for d in parents for p in sorted((ROOT.parent/d).rglob('*')) if p.is_file() and '__pycache__' not in p.parts};js('authority/parent_authority.json',{'files':files,'snapshot_utc':utc()})
    manifest=json.loads((P4/'completion/artifact_manifest.json').read_text())['files']
    for p,r in manifest.items():assert sha(P4/p)==r['sha256']
    receipt=json.loads((P4/'reviewer_release/verification/validation_receipt.json').read_text());assert receipt['passed']
    js('authority/parent_validation.json',{'passed':True,'parent_manifest_verified':len(manifest),'parent_full_replay_receipt':receipt,'full_rebuild_performed':False,'parent_immutable_files':len(files)})
    shutil.copyfile(directive,ROOT/'authority/research_directive.txt');part=read(P4/'authority/feature_partition.csv.gz');gs=[{k:r[k] for k in ('geneId','geneIndex')} for r in part if r['panel']=='state'];gr=[{k:r[k] for k in ('geneId','geneIndex')} for r in part if r['panel']=='realization']
    assert len(gs)==14660 and len(gr)==6339 and not set(r['geneId'] for r in gs)&set(r['geneId'] for r in gr)
    write('authority/primary_GS.csv.gz',gs);write('authority/primary_GR.csv.gz',gr);shutil.copyfile(P4/'authority/feature_partition.csv.gz',ROOT/'authority/feature_partition.csv.gz')
    js('authority/feature_partition_receipt.json',{'source':'relational_morphogenesis12/genes/state_realization_partition.csv.gz','partition_sha256':sha(ROOT/'authority/feature_partition.csv.gz'),'G_S_sha256':sha(ROOT/'authority/primary_GS.csv.gz'),'G_R_sha256':sha(ROOT/'authority/primary_GR.csv.gz'),'partition_seed':'RMMO12_STATE_REALIZATION_V1','gene_order_preserved':True,'disjoint':True})
    raw=sp.load_npz(P4/'normalization/selected_eligible_counts.npz');idx=np.load(P4/'normalization/count_index.npz');state=raw[:,idx['state_columns']].copy();np.testing.assert_array_equal(idx['geneId'][idx['state_columns']],[r['geneId'] for r in gs]);(ROOT/'data_authority').mkdir(exist_ok=True);sp.save_npz(ROOT/'data_authority/GS_raw_counts.npz',state)
    cells=read(P4/'authority/cell_carrier.csv.gz');write('data_authority/cell_identity.csv.gz',[{'cellHash':r['cellHash'],'replicate':r['replicate']} for r in cells]);assert len(cells)==3704
    js('data_authority/preparation_receipt.json',{'phase':'preparation OUTSIDE Stage A','source_selected_counts_sha256':sha(P4/'normalization/selected_eligible_counts.npz'),'GS_count_input_sha256':sha(ROOT/'data_authority/GS_raw_counts.npz'),'GS_shape':list(state.shape),'source_all_gene_carrier_opened_here':True,'GR_or_annotations_supplied_to_StageA':False,'new_strict_outcomes_inspected':False})
    write('authority/construction_dependency_registry.csv',[{'variable':v,'source':s,'reason':r,'permitted':True} for v,s,r in [('GS_counts','data_authority/GS_raw_counts.npz','only ordered frozen G_S columns'),('GS_geneId_geneIndex','authority/primary_GS.csv.gz','membership and signed-hash projection'),('cellHash_replicate','data_authority/cell_identity.csv.gz','ordered cell alignment and pair carriers'),('plan_seeds_constants','freeze/analysis_plan.json and strict_lane_plan.csv','precommitted deterministic lane/geometry specification')]])
    write('authority/forbidden_input_registry.csv',[{'input':v,'construction_permitted':False,'enforcement':'file allowlist, independent membership audit, poison replay and no evaluation imports'} for v in freeze()['forbidden_construction_inputs']])
    if manuscript:
        out=ROOT/'manuscript_authority';out.mkdir(exist_ok=True);shutil.copyfile(manuscript,out/'Paper3_v4_original.tex');js('manuscript_authority/source_receipt.json',{'source_filename':Path(manuscript).name,'sha256':sha(out/'Paper3_v4_original.tex'),'original_modified':False,'new_v5_output_scoped_to_this_goal':True})
    print('Preparation complete; no strict outcomes evaluated; parents',len(files),flush=True)
def audit_genes():
    # Independent auditor may inspect G_R identity, but never supplies it to construction.
    gs=read(ROOT/'authority/primary_GS.csv.gz');gr={r['geneId'] for r in read(ROOT/'authority/primary_GR.csv.gz')};s={r['geneId'] for r in gs};audit=[];union=set();membership=[]
    for spec in read(ROOT/'freeze/strict_lane_plan.csv'):
        genes=[r for r in gs if not spec['panel_seed'] or stable(spec['panel_seed']+'|'+r['geneId'])%10<7];g={r['geneId'] for r in genes};assert len(genes)==int(spec['gene_count']) and g<=s and not g&gr;union|=g
        audit.append({'lane':spec['lane'],'gene_count':len(g),'intersection_with_GS':len(g&s),'intersection_with_GR':len(g&gr),'ordered_membership_sha256':hashlib.sha256('|'.join(r['geneId'] for r in genes).encode()).hexdigest()})
        membership.extend({'lane':spec['lane'],'lane_gene_position':i,**r} for i,r in enumerate(genes))
    assert union<=s and not union&gr;write('strict_lanes/lane_gene_membership.csv.gz',membership);write('strict_lanes/lane_intersection_audit.csv',audit)
    js('separation_audit/union_gene_audit.json',{'union_genes':len(union),'union_subset_GS':True,'intersection_GR':0,'all_nine_verified':True,'auditor_not_construction':True,'audit_utc':utc()});print('Nine strict gene sets certified; union = G_S, GR intersection 0',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['seal','prepare','audit_genes']);p.add_argument('--directive');p.add_argument('--manuscript');a=p.parse_args();{'seal':seal,'prepare':lambda:prepare(a.directive,a.manuscript),'audit_genes':audit_genes}[a.action]()

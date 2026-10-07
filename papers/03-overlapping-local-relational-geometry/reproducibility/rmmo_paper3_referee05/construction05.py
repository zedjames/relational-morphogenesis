"""Stage A: completely G_S-only numerical construction, no evaluation imports."""
import builtins,io,contextlib
from collections import Counter
import scipy.sparse as sp
from io05 import *
INPUTS=('data_authority/GS_raw_counts.npz','authority/primary_GS.csv.gz','data_authority/cell_identity.csv.gz','freeze/analysis_plan.json','freeze/strict_lane_plan.csv')
@contextlib.contextmanager
def access_guard():
    allowed={str((ROOT/p).resolve()) for p in INPUTS};events=[];orig,iorig=builtins.open,io.open
    def checked(original):
        def op(file,*args,**kwargs):
            if isinstance(file,(str,Path)):
                path=str(Path(file).resolve());mode=args[0] if args else kwargs.get('mode','r');assert path in allowed and not any(x in mode for x in 'wa+'),('Stage A forbidden access',path,mode);events.append({'path':str(Path(path).relative_to(ROOT)),'mode':mode})
            return original(file,*args,**kwargs)
        return op
    builtins.open,io.open=checked(orig),checked(iorig)
    try:yield events
    finally:builtins.open,io.open=orig,iorig
def compute(raw,genes,cells,plan,specs):
    assert raw.shape==(3704,14660) and set(cells[0])=={'cellHash','replicate'} and set(genes[0])=={'geneId','geneIndex'}
    assert np.isfinite(raw.data).all() and np.all(raw.data>=0);ids=[g['geneId'] for g in genes];ix=[np.flatnonzero([c['replicate']==r for c in cells]) for r in REPS];arrays={};lanes={};registry=[];thresholds=[];votes={p:Counter() for p in PAIRS}
    for spec in specs:
        lane=spec['lane'];cols=np.asarray([i for i,g in enumerate(ids) if not spec['panel_seed'] or stable(spec['panel_seed']+'|'+g)%10<7]);assert len(cols)==int(spec['gene_count'])
        x=raw[:,cols];lib=np.asarray(x.sum(1)).ravel();x=sp.diags(np.divide(10000.,lib,out=np.zeros_like(lib,dtype=float),where=lib>0))@x.astype(float);x=x.tocsr();x.data=np.log1p(x.data)
        buckets=np.array([stable(spec['projection_seed']+'|bucket|'+ids[i])%64 for i in cols]);signs=np.array([1 if stable(spec['projection_seed']+'|sign|'+ids[i])%2==0 else -1 for i in cols]);h=sp.csc_matrix((signs,(np.arange(len(cols)),buckets)),shape=(len(cols),64));v=unit((x@h).toarray());v=unit(v.astype(np.float32)) if lane=='S0' else v
        assert np.isfinite(v).all();lanes[lane]=v;arrays['lane_'+lane+'_embedding']=v;arrays['lane_'+lane+'_library']=lib;arrays['lane_'+lane+'_columns']=cols
        registry.append({'lane':lane,'role':spec['role'],'genes':len(cols),'projection_seed':spec['projection_seed'],'panel_seed':spec['panel_seed'],'zero_library_cells':int(np.sum(lib==0)),'own_denominator_only':True,'embedding_digest':array_digest(v)})
        for pair,(a,b) in ENDS.items():
            d=chord(v[ix[a]],v[ix[b]]);eps=float(np.quantile(np.r_[d.min(1),d.min(0)],.5));li,ti=np.where(d<=eps);left,right=ix[a][li],ix[b][ti];votes[pair].update((left*3704+right).tolist())
            arrays[lane+'_'+pair+'_left']=left;arrays[lane+'_'+pair+'_right']=right;arrays[lane+'_'+pair+'_epsilon']=np.asarray(eps);thresholds.append({'lane':lane,'pair':pair,'epsilon':eps,'edges':len(left)})
    z=lanes['S0'];summaries=[];edges=[];geo={};core={}
    for pair,(a,b) in ENDS.items():
        codes=np.array(sorted(c for c,n in votes[pair].items() if n>=7),dtype=np.int64);left,right=codes//3704,codes%3704;assert len(left)>0,('empty strict pair: do not manufacture geometry',pair);core[pair]=(left,right);arrays[pair+'_left']=left;arrays[pair+'_right']=right
        summaries.append({'pair':pair,'edges':len(left),'source_endpoints':len(set(left)),'target_endpoints':len(set(right)),'source_max_degree':max(Counter(left).values()),'target_max_degree':max(Counter(right).values())});edges.extend({'pair':pair,'source':int(l),'target':int(t),'lane_support':votes[pair][int(l*3704+t)]} for l,t in zip(left,right))
        for a0,b0,l,t in [(a,b,left,right),(b,a,right,left)]:
            name=f'L{a0+1}{b0+1}';anchors=np.unique(l);targets=np.unique(t);d=chord(z[ix[a0]],z[anchors]);rho=d.min(1);adj=np.zeros((len(anchors),len(targets)),dtype=np.int32);adj[np.searchsorted(anchors,l),np.searchsorted(targets,t)]=1;geo[name]={};arrays[name+'_anchors']=anchors
            for q,quant in SCALES.items():
                rad=float(np.quantile(rho,quant));nodes=ix[a0][rho<=rad];active=d[rho<=rad]<=rad;mem=active.astype(np.int32)@adj>0;assert mem.any(1).all()
                arrays[name+'_'+q+'_nodes']=nodes;arrays[name+'_'+q+'_radius']=np.asarray(rad);arrays[name+'_'+q+'_active']=active;arrays[name+'_'+q+'_fiber_targets']=targets;arrays[name+'_'+q+'_fiber_membership']=mem;geo[name][q]=nodes
    for name,a,b in OVERLAPS:
        for q in SCALES:arrays[name+'_'+q]=np.intersect1d(geo[a][q],geo[b][q])
    return arrays,registry,thresholds,summaries,edges
def load_and_compute():
    # Import/decompression modules are preloaded outside the enforced construction boundary.
    with access_guard() as events:
        raw=sp.load_npz(ROOT/INPUTS[0]);genes=read(ROOT/INPUTS[1]);cells=read(ROOT/INPUTS[2]);plan=json.loads((ROOT/INPUTS[3]).read_text());specs=read(ROOT/INPUTS[4]);result=compute(raw,genes,cells,plan,specs)
    return result,events
def main():
    freeze();audit=json.loads((ROOT/'separation_audit/union_gene_audit.json').read_text());assert audit['intersection_GR']==0 and audit['all_nine_verified']
    result,events=load_and_compute();a,reg,th,summ,edges=result;save('geometry/stageA_arrays.npz',**a);save('geometry/pair_cores.npz',**{p+'_'+s:a[p+'_'+s] for p in PAIRS for s in ('left','right')})
    for lane in ('S'+str(i) for i in range(9)):
        save(f'strict_lanes/lane_embeddings/{lane}.npz',embedding=a['lane_'+lane+'_embedding'],library=a['lane_'+lane+'_library'],GS_columns=a['lane_'+lane+'_columns']);save(f'strict_lanes/lane_relations/{lane}.npz',**{k:v for k,v in a.items() if k.startswith(lane+'_')})
    write('strict_lanes/lane_registry.csv',reg);write('strict_lanes/lane_thresholds.csv',th);write('geometry/pair_landmark_summary.csv',summ);write('geometry/pair_relations.csv.gz',edges);write('strict_lanes/consensus_edges.csv.gz',edges)
    for q in SCALES:write(f'geometry/{q}_domains.csv.gz',[{'orientation':n,'cell':int(v),'radius':float(a[n+'_'+q+'_radius'])} for n in ('L12','L21','L13','L31','L23','L32') for v in a[n+'_'+q+'_nodes']])
    write('geometry/overlaps.csv',[{'overlap':n,'scale':q,'cells':len(a[n+'_'+q])} for n,_,_ in OVERLAPS for q in SCALES])
    hashes={k:array_digest(v) for k,v in a.items()};js('separation_audit/geometry_pre_reveal_hashes.json',{'array_content_hashes':hashes,'container_sha256':sha(ROOT/'geometry/stageA_arrays.npz'),'GS_only':True,'geometry_utc':utc(),'heldout_loaded':False})
    js('separation_audit/construction_input_manifest.json',{'files':{p:sha(ROOT/p) for p in INPUTS},'forbidden_file_access_raises':True,'construction_module_sha256':sha(Path(__file__))})
    js('separation_audit/heldout_access_audit.json',{'GR_counts_read_during_StageA':False,'ZR_loaded_during_StageA':False,'GR_columns_used':0,'annotations_used':False,'history_used':False,'heldout_statistics_accessed':False,'instrumented_file_accesses':events,'enforcement':'both builtins.open and io.open allowlisted; no parent/evaluation imports; pure computation receives only GS arrays and identity'})
    js('strict_lanes/consensus_receipt.json',{'lanes':9,'criterion':7,'calibration':'pairwise q50 quantile of all directional nearest distances','union_GR_intersection':0,'consensus_sha256':sha(ROOT/'strict_lanes/consensus_edges.csv.gz')});js('geometry/geometry_receipt.json',{'state_only':True,'fibers_nonempty':'CONSTRUCTION_GUARANTEE','array_digest_count':len(hashes),'primary':'strict GS-only geometry','finished_utc':utc()})
    print('Strict Stage A complete',summ,[(n,len(a[n+'_q25'])) for n,_,_ in OVERLAPS],flush=True)
if __name__=='__main__':main()

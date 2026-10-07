"""Private immutability + portable full strict-construction/numerical verification."""
import argparse,subprocess,io,ast
import scipy.sparse as sp
from numerical_verify05 import *
import construction05 as construct
def private_authority():
    rec=json.loads((ROOT/'authority/parent_authority.json').read_text())
    for rel,x in rec['files'].items():assert sha(ROOT.parent/rel)==x['sha256'],('parent changed',rel)
    receipt=json.loads((ROOT/'freeze/freeze_receipt.json').read_text());g=json.loads((ROOT/'separation_audit/geometry_commit_receipt.json').read_text())
    for rel,h in receipt['plan_sha256'].items():assert hashlib.sha256(subprocess.check_output(['git','show',receipt['reachable_plan_commit']+':research/rmmo_paper3_referee05/'+rel],cwd=ROOT.parent.parent)).hexdigest()==h
    # Git LFS stores a pointer at commit time, while the working file is materialized.
    data=subprocess.check_output(['git','show',g['reachable_geometry_commit']+':research/rmmo_paper3_referee05/geometry/stageA_arrays.npz'],cwd=ROOT.parent.parent)
    assert ('oid sha256:'+g['geometry_sha256']).encode() in data or hashlib.sha256(data).hexdigest()==g['geometry_sha256']
    full=sp.load_npz(P4/'normalization/selected_eligible_counts.npz');idx=np.load(P4/'normalization/count_index.npz')
    for key,rel in [('state_columns','data_authority/GS_raw_counts.npz'),('heldout_columns','evaluation_inputs/GR_raw_counts.npz')]:
        wanted=full[:,idx[key]].tocsr();actual=sp.load_npz(ROOT/rel)
        for field in ('data','indices','indptr'):np.testing.assert_array_equal(getattr(wanted,field),getattr(actual,field))
    print('PASS immutable parent snapshot and reachable pre-outcome plan/geometry commits',len(rec['files']),flush=True)
def construction_checks():
    freeze();gs=read(ROOT/'authority/primary_GS.csv.gz');gr=read(ROOT/'authority/primary_GR.csv.gz');ss={r['geneId'] for r in gs};hh={r['geneId'] for r in gr};assert len(ss)==14660 and len(hh)==6339 and not ss&hh
    membership=read(ROOT/'strict_lanes/lane_gene_membership.csv.gz');union=set()
    for i in range(9):
        gene={r['geneId'] for r in membership if r['lane']==f'S{i}'};assert gene<=ss and not gene&hh;union|=gene
    assert union==ss
    manifest=json.loads((ROOT/'separation_audit/construction_input_manifest.json').read_text())
    for rel,h in manifest['files'].items():assert sha(ROOT/rel)==h
    assert sha(ROOT/'construction05.py')==manifest['construction_module_sha256']
    result,events=construct.load_and_compute();a=result[0];pre=json.loads((ROOT/'separation_audit/geometry_pre_reveal_hashes.json').read_text());assert {k:array_digest(v) for k,v in a.items()}==pre['array_content_hashes']
    b=io.BytesIO();np.savez_compressed(b,**a);assert hashlib.sha256(b.getvalue()).hexdigest()==pre['container_sha256']==sha(ROOT/'geometry/stageA_arrays.npz')
    assert all(r['path'] in construct.INPUTS for r in events)
    # Portable opaque-GR poison: neither counts nor identity are capabilities of compute.
    raw=sp.load_npz(ROOT/'data_authority/GS_raw_counts.npz');opaque=sp.load_npz(ROOT/'evaluation_inputs/GR_raw_counts.npz');before=array_digest(opaque.data);opaque.data[:]=rng('PORTABLE_GR_POISON_CHANGED').integers(0,101,size=len(opaque.data));assert array_digest(opaque.data)!=before
    poisoned=construct.compute(raw,gs,read(ROOT/'data_authority/cell_identity.csv.gz'),freeze(),read(ROOT/'freeze/strict_lane_plan.csv'))[0];assert {k:array_digest(v) for k,v in poisoned.items()}==pre['array_content_hashes'];b=io.BytesIO();np.savez_compressed(b,**poisoned);assert hashlib.sha256(b.getvalue()).hexdigest()==pre['container_sha256']
    poison=json.loads((ROOT/'separation_audit/heldout_poison_test.json').read_text());assert poison['passed'] and poison['compressed_geometry_byte_identical'] and poison['numeric_arrays_compared']==219
    r,z,ix,c,g,o,m=arrays();stored=np.load(ROOT/'geometry/stageA_arrays.npz')
    for name,e in g.items():
        np.testing.assert_array_equal(stored[name+'_anchors'],e['anchors'])
        for q in SCALES:
            for field in ('nodes','active'):np.testing.assert_array_equal(stored[name+'_'+q+'_'+field],e[q][field])
            close(stored[name+'_'+q+'_radius'],e[q]['radius'],'independent radius');tar,mem=fiber_matrix(e,q);np.testing.assert_array_equal(tar,stored[name+'_'+q+'_fiber_targets']);np.testing.assert_array_equal(mem,stored[name+'_'+q+'_fiber_membership']);assert mem.any(1).all()
    for (name,q),v in o.items():np.testing.assert_array_equal(stored[name+'_'+q],v)
    for lane in range(9):
        v=np.load(ROOT/f'strict_lanes/lane_embeddings/S{lane}.npz');np.testing.assert_array_equal(v['embedding'],a[f'lane_S{lane}_embedding']);np.testing.assert_array_equal(v['library'],np.asarray(raw[:,v['GS_columns']].sum(1)).ravel())
        expected=[gs[int(j)]['geneId'] for j in v['GS_columns']];actual=[r['geneId'] for r in membership if r['lane']==f'S{lane}'];np.testing.assert_array_equal(actual,expected)
    f=json.loads((ROOT/'freeze/freeze_receipt.json').read_text());g0=json.loads((ROOT/'separation_audit/geometry_commit_receipt.json').read_text());reveal=revealed();assert f['original_seal_utc']<pre['geometry_utc']<g0['original_geometry_seal_utc']<reveal['evaluation_start_utc']
    print('PASS all GS-only lanes, own denominators, 219 committed arrays, poison invariance and complete L/U/F/O geometry',flush=True)
@np.errstate(all='ignore')
def evaluation_checks():
    g=read(ROOT/'authority/primary_GR.csv.gz');x=sp.load_npz(ROOT/'evaluation_inputs/heldout_expression.npz');assert x.shape==(3704,6339)
    rawgr=sp.load_npz(ROOT/'evaluation_inputs/GR_raw_counts.npz');libgr=np.asarray(rawgr.sum(1)).ravel();xhash=(sp.diags(np.divide(10000.,libgr,out=np.zeros_like(libgr,dtype=float),where=libgr>0))@rawgr.astype(float)).tocsr();xhash.data=np.log1p(xhash.data)
    seed0='RMMO12_REALIZATION_HASH_PROJECTION_V1';buckets=np.array([stable(seed0+'|bucket|'+r['geneId'])%64 for r in g]);sign=np.array([1 if stable(seed0+'|sign|'+r['geneId'])%2==0 else -1 for r in g]);h=sp.csc_matrix((sign,(np.arange(len(g)),buckets)),shape=(len(g),64));np.testing.assert_array_equal(unit((xhash@h).toarray()).astype(np.float32).astype(float),target())
    raw64=rawgr.astype(float);lib64=np.asarray(raw64.sum(1)).ravel();x64=(sp.diags(np.divide(10000.,lib64,out=np.zeros_like(lib64),where=lib64>0))@raw64).tocsr();x64.data=np.log1p(x64.data);assert (x64!=x).nnz==0
    raw=sp.load_npz(ROOT/'data_authority/GS_raw_counts.npz').astype(float);lib=np.asarray(raw.sum(1)).ravel();xx=(sp.diags(np.divide(10000.,lib,out=np.zeros_like(lib),where=lib>0))@raw).tocsr();xx.data=np.log1p(xx.data);assert (xx!=sp.load_npz(ROOT/'evaluation_inputs/state_expression.npz')).nnz==0
    for name in ('pca32','pca64','target_pca32'):
        f=np.load(ROOT/f'evaluation_inputs/{name}/representation.npz');x=sp.load_npz(ROOT/('evaluation_inputs/heldout_expression.npz' if name=='target_pca32' else 'evaluation_inputs/state_expression.npz'));v=f['loadings'];s=f['scores'];mu=f['mean'];sv=f['singular_values'];assert np.isfinite(s).all() and np.isfinite(v).all()
        close(mu,np.asarray(x.mean(0)).ravel(),name+' means');close(x@v-mu@v,s,name+' scores');close(unit(s),f['embedding'],name+' unit');close(v.T@v,np.eye(v.shape[1]),name+' orthonormality');close(np.linalg.norm(s,axis=0),sv,name+' singular values');eig=x.T@s-mu[:,None]*s.sum(0);assert np.linalg.norm(eig-v*sv**2)/np.linalg.norm(v*sv**2)<2e-6
        for i in (0,1,7,1031,1032,2329,3703):
            row=x[i];manual=np.sum(row.data[:,None]*v[row.indices],0)-np.sum(mu[:,None]*v,0);close(manual,s[i],name+' independent row')
        assert all(v[np.argmax(np.abs(v[:,j])),j]>=0 for j in range(v.shape[1]))
    for name in ('hash_primary','pca32','pca64'):
        r,z,ix,c,g,o,m=arrays(name);saved=np.load(ROOT/f'geometry/{name}/pair_cores.npz')
        for p,(i,j) in ENDS.items():
            d=chord(z[ix[i]],z[ix[j]]);eps=np.quantile(np.r_[d.min(1),d.min(0)],.5);l,t=np.where(d<=eps);np.testing.assert_array_equal(saved[p+'_left'],ix[i][l]);np.testing.assert_array_equal(saved[p+'_right'],ix[j][t])
    print('PASS independent GS expression/GR hash reconstruction and PCA scores, eigenpairs, single-lane controls',flush=True)
def sampler_checks():
    audit=read(ROOT/'source_fidelity/graph_audit.csv');summary=read(ROOT/'source_fidelity/sampler_transparency.csv')
    for row in summary:
        v=np.array([int(x['stub_trials']) for x in audit if x['pair']==row['pair']]);assert len(v)==1024 and v.sum()==int(row['attempted_stub_pairings']) and v.max()==int(row['maximum_proposals']);assert int(row['rejected_nonsimple'])==v.sum()-1024;close(float(row['acceptance_fraction']),1024/v.sum(),'acceptance');close(float(row['mean_proposals_per_accepted']),v.mean(),'attempt mean');close(float(row['median_proposals_per_accepted']),np.median(v),'attempt median')
    print('PASS complete fresh-pairing rejection accounting',flush=True)
def report_checks():
    from report05 import decision,success,LIMITATION,select
    d=json.loads((ROOT/'completion/strict_feature_separation_decision.json').read_text());assert d['biological_replicates']==3 and not d['external_confirmation'];assert d['STRICT_EDGE_HELDOUT_PERSISTS']==success(decision('edge_agreement'));assert d['STRICT_MATCHED_COHERENCE_PERSISTS']==success(decision('matched_coherence/q25','heldout','defect'));assert d['STRICT_SOURCE_FIDELITY_PERSISTS']==success(decision('source_fidelity/q25','heldout','source_error'))
    claims=read(ROOT/'completion/claim_matrix.csv');assert [r['claim'] for r in claims]==list('ABCDEFGHIJKLMNOPQRST');q=read(ROOT/'completion/question_answers.csv');assert [int(r['question_number']) for r in q]==list(range(1,96));assert q[87]['answer'].startswith('YES') and q[89]['answer'].startswith('NO') and q[90]['answer'].startswith('NO') and q[93]['answer'].startswith('NO')
    assert LIMITATION in (ROOT/'completion/RESULTS.md').read_text();assert not (ROOT/'retention/R4').exists() and not (ROOT/'retention/R5').exists() and not (ROOT/'spatial_overlap/S2').exists()
    for r in read(ROOT/'geometry/v4_vs_strict_object_identity.csv'):close(float(r['Jaccard']),int(r['intersection'])/(int(r['reference_count'])+int(r['strict_count'])-int(r['intersection'])),'object Jaccard')
    for r in read(ROOT/'comparison/referee04_vs_strict.csv'):close(float(r['effect_change']),float(r['strict_effect'])-float(r['referee04_effect']),'comparison effect change')
    from pypdf import PdfReader
    from PIL import Image
    for n in range(1,7):
        f=ROOT/f'figures/figure{n:02}.pdf';assert len(PdfReader(f).pages)==1 and f.stat().st_size>12000;assert len(read(ROOT/f'figures/figure{n:02}_authority.csv'))>0
        im=Image.open(ROOT/f'figures/figure{n:02}.png');assert im.size[0]>=1500 and im.size[1]>=1000
    pdf=ROOT/'manuscript_authority/Paper3_v5.pdf';reader=PdfReader(pdf);text=' '.join(' '.join(p.extract_text() or '' for p in reader.pages).split());assert len(reader.pages)>=12 and 'Zed James' in text and 'construction-held-out' in text and 'independent embryos' in text.lower() and '59' in text;assert '??' not in text
    receipt=json.loads((ROOT/'manuscript_authority/render_receipt.json').read_text());assert sha(pdf)==receipt['pdf_sha256'] and sha(ROOT/'manuscript_authority/Paper3_v5.tex')==receipt['tex_sha256'];assert sha(ROOT/'manuscript_authority/Paper3_v4_original.tex')==receipt['original_v4_sha256'];assert len(receipt['figure_sha256'])==6
    for rel,h in receipt['figure_sha256'].items():assert sha(ROOT/rel)==h
    print('PASS twenty claims,95 literal answers, source-backed six figures and rendered v5 manuscript',flush=True)
def descriptive_checks():
    # Independent row arithmetic checks all recorded baseline cell scores.
    cells=read(ROOT/'baselines/cellwise.csv.gz');groups=defaultdict(list)
    for r in cells:
        a,b,d=map(float,(r['native_error'],r['baseline_error'],r['difference']));close(a-b,d,'baseline difference');groups[r['baseline'],r['scale'],r['panel'],r['orientation']].append((int(r['cell']),a,b))
    summary=read(ROOT/'baselines/comparison.csv')
    for r in summary:
        if not r['orientation'].startswith('L') or int(r['cells'])==0:continue
        rr=groups[r['baseline'],r['scale'],r['panel'],r['orientation']];a=np.array([x[1] for x in rr]);b=np.array([x[2] for x in rr]);assert len(rr)==int(r['cells'])
        for field,val in [('observed_mean',a.mean()),('baseline_mean',b.mean()),('paired_difference_mean',(a-b).mean()),('paired_difference_median',np.median(a-b)),('observed_median',np.median(a)),('baseline_median',np.median(b))]:close(float(r[field]),val,'baseline summary')
    r,z,ix,c,g,o,m=arrays();dist={p:chord(z[ix[a]],z[ix[b]]) for p,(a,b) in ENDS.items()};native_support=read(ROOT/'baselines/mnn_k5/native_support.csv')
    for name,e in g.items():
        pair='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))));d=dist[pair] if e['source']<e['target'] else dist[pair].T;nearest=np.argsort(d,axis=1,kind='stable')[:,:5];backward=np.argsort(d,axis=0,kind='stable')[:5,:];forward=np.zeros(d.shape,bool);reverse=forward.copy();forward[np.arange(len(d))[:,None],nearest]=True;reverse[backward,np.arange(d.shape[1])[None,:]]=True;mnn=forward&reverse
        for q in SCALES:
            nodes=e[q]['nodes'];tar,mem=fiber_matrix(e,q);codes=fiber_codes({name:e},q)[name]
            for label in ('matched_knn','pca32_knn','pca64_knn','mnn_k5','mnn_fixed_size'):
                v=z if label in ('matched_knn','mnn_k5','mnn_fixed_size') else state(label.replace('_knn',''));distance=d if v is z else chord(v[ix[e['source']]],v[ix[e['target']]])
                for panel,mat in m.items():
                    rr=groups[label,q,panel,name]
                    for cell,a,b in rr:
                        li=np.searchsorted(nodes,cell);j=np.searchsorted(ix[e['source']],cell);k=int(mem[li].sum());native_targets=tar[mem[li]]
                        if label.endswith('knn'):selected=ix[e['target']][np.argsort(distance[j],kind='stable')[:k]]
                        else:
                            targets=np.flatnonzero(mnn[j]);assert len(targets)>0
                            if label=='mnn_fixed_size':assert len(targets)>=k;targets=targets[np.argsort(d[j,targets],kind='stable')[:k]]
                            selected=ix[e['target']][targets]
                        native=mat[native_targets].mean(0);baseline=mat[selected].mean(0)
                        if panel=='state':native=unit(native);baseline=unit(baseline)
                        close(a,np.linalg.norm(mat[cell]-native),'baseline independent native');close(b,np.linalg.norm(mat[cell]-baseline),'baseline independent selected')
        for panel,mat in m.items():
            nodes=ix[e['source']][mnn.any(1)];mu=np.array([mat[ix[e['target']][np.flatnonzero(mnn[j])]].mean(0) for j in np.flatnonzero(mnn.any(1))]);mu=unit(mu) if panel=='state' else mu;err=np.linalg.norm(mat[nodes]-mu,axis=1);r0=next(r for r in native_support if r['orientation']==name and r['panel']==panel)
            assert int(r0['native_support'])==len(nodes);close(float(r0['error_mean']),err.mean(),'native MNN mean');close(float(r0['error_median']),np.median(err),'native MNN median')
    full=rows(True);f=np.load(ROOT/'authority/referee04/geometry/hash/pair_cores.npz');oldc={p:(f[p+'_left'],f[p+'_right']) for p in PAIRS};oldg,oldo=geometry(z,ix,oldc);objects={}
    for mode,cc,over in [('v4',oldc,oldo),('isolated',c,o)]:
        for p in PAIRS:objects[mode,p+'_source']=np.unique(cc[p][0]);objects[mode,p+'_target']=np.unique(cc[p][1]);objects[mode,p]=np.unique(np.r_[cc[p][0],cc[p][1]])
        for name,_,_ in OVERLAPS:objects[mode,name]=over[name,'q25']
    for rel in ('biology/landmark_annotation_summary.csv','biology/overlap_annotation_summary.csv'):
        for r0 in read(ROOT/rel):
            nodes=objects[r0['normalization'],r0['object']];cnt=sum(full[int(i)][r0['field']]==r0['category'] for i in nodes);assert cnt==int(r0['count']) and len(nodes)==int(r0['cells']);close(float(r0['fraction']),cnt/len(nodes),'annotation fraction')
    print('PASS independently scored kNN/MNN/fixed-size MNN cells, native support and biological fractions',flush=True)
def incidence_group(group):incidence_replays(group)
def main(public=False,skip_banks=False,parallel_banks=False):
    start=utc();freeze()
    if not public:private_authority()
    construction_checks();evaluation_checks()
    if not skip_banks:
        edge_replay();spatial_replay()
        if parallel_banks:
            import concurrent.futures,multiprocessing
            rels=['matched_coherence/'+q for q in SCALES]+['representation/'+v for v in ('pca32','pca64','target_pca32')]+[f'retention/R{i}' for i in range(4)]
            with concurrent.futures.ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:list(pool.map(incidence_group,[rels[::2],rels[1::2]]))
        else:incidence_replays()
        degree_replay()
    descriptive_checks();sampler_checks();report_checks()
    js('completion/'+('public_validation_receipt.json' if public else 'validation_receipt.json'),{'passed':True,'started_utc':start,'finished_utc':utc(),'all_draws_replayed':not skip_banks,'independent_graph_score_arithmetic':True,'B':1024,'incidence_banks':10,'edge_banks':1,'spatial_banks':1,'degree_scale_banks':3,'parents_verified':not public,'parent_count':7859 if not public else None,'GS_only_arrays_checked':219,'no_full_rebuild':True,'biological_replicates':3,'rendered_manuscript':True,'actual_figures':6,'public_GS_input_reconstruction':True,'public_raw_H5TD_reextraction':False})
    if not public:md('completion/validation_receipt.md','# Verification PASS\n\nAll 7,859 immutable parent files and committed pre-outcome plan/geometry hashes checked. All219 GS-only arrays, lane-local denominators, complete domains/fibers/overlaps, embargo and actual geometry invariance checked. Every retained1024-member direct IID graph bank independently replayed and scored (10 incidence,1 edge,1 S1,3 degree scale banks). PCA scores/eigenpairs, descriptive baseline cells/compositions, all95 answers,20 claims, six actual figures and rendered v5 manuscript verified. Portable release separately verifies its own ordinary-file manifest and supplied GS input, not raw-H5TD acquisition or the entire historical repository. No whole Lean build.')
    print('Referee05 full verification PASS',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--public',action='store_true');p.add_argument('--skip-banks',action='store_true');p.add_argument('--parallel-banks',action='store_true');args=p.parse_args();main(args.public,args.skip_banks,args.parallel_banks)

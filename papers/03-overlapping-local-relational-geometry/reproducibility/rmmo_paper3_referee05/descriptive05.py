"""Descriptive baseline/representation/annotation updates, no population tests."""
import argparse
from scipy.stats import spearmanr
from common05 import *
from objects import errors,summary_difference,correlation
from controls05 import compare
def baselines():
    freeze();r,z,ix,c,g,o,m=arrays();results=[];cells=[];support=[];embs={'matched_knn':z,'pca32_knn':state('pca32'),'pca64_knn':state('pca64')}
    distances={name:{p:chord(emb[ix[a]],emb[ix[b]]) for p,(a,b) in ENDS.items()} for name,emb in embs.items()}
    mnn={};native_mnn=[]
    for name,e in g.items():
        p='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))));d=distances['matched_knn'][p];d=d if e['source']<e['target'] else d.T
        nearest=np.argsort(d,axis=1,kind='stable')[:,:5];reverse=np.argsort(d,axis=0,kind='stable')[:5,:]
        f=np.zeros(d.shape,bool);f[np.arange(len(d))[:,None],nearest]=True;back=np.zeros(d.shape,bool);back[reverse,np.arange(d.shape[1])[None,:]]=True;member=f&back
        mnn[name]=(d,member)
        # Native MNN support is independent of pair-native source domains.
        # Do not substitute common-support comparisons for this separate estimand.
        available=ix[e['source']][member.any(1)]
        selected=[ix[e['target']][np.flatnonzero(member[j])] for j in np.flatnonzero(member.any(1))]
        for panel,mat in m.items():
            mu=np.asarray([mat[t].mean(0) for t in selected]);mu=unit(mu) if panel=='state' else mu
            err=np.linalg.norm(mat[available]-mu,axis=1)
            native_mnn.append({'baseline':'mnn_k5','orientation':name,'panel':panel,'native_support':len(available),'error_mean':float(err.mean()),'error_median':float(np.median(err)),
                'mean_target_size':float(np.mean([len(t) for t in selected])),'descriptive_only':True,'pair_native_comparison':'not defined outside common support'})
    for q in SCALES:
        native=observed_means(g,q,m);ne,nd=errors(native,g,o,q,m)
        for label in tuple(embs)+('mnn_k5','mnn_fixed_size'):
            by={}
            for name,e in g.items():
                nodes=e[q]['nodes'];tar,mem=fiber_matrix(e,q);sizes=mem.sum(1);p='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))))
                if label in embs:
                    d=distances[label][p];d=d if e['source']<e['target'] else d.T;order=np.argsort(d[np.searchsorted(ix[e['source']],nodes)],axis=1,kind='stable')
                    selected=[ix[e['target']][order[i,:k]] for i,k in enumerate(sizes)];take=np.arange(len(nodes));native_support=len(nodes)
                else:
                    d,member=mnn[name];rows_local=np.searchsorted(ix[e['source']],nodes);take=[];selected=[]
                    for i,j in enumerate(rows_local):
                        targets=np.flatnonzero(member[j])
                        if not len(targets) or label=='mnn_fixed_size' and len(targets)<sizes[i]:continue
                        if label=='mnn_fixed_size':targets=targets[np.argsort(d[j,targets],kind='stable')[:sizes[i]]]
                        take.append(i);selected.append(ix[e['target']][targets])
                    take=np.asarray(take,int);native_support=int(member.any(1).sum())
                available=nodes[take];support.append({'baseline':label,'scale':q,'orientation':name,'native_baseline_support':native_support,'pair_native_support':len(nodes),'common_support':len(available),
                    'size_matched':label!='mnn_k5','mean_target_size':float(np.mean([len(t) for t in selected])) if selected else 0.})
                for panel,mat in m.items():
                    mu=np.asarray([mat[t].mean(0) for t in selected]);mu=unit(mu) if panel=='state' and len(mu) else mu
                    ev=np.linalg.norm(mat[available]-mu,axis=1) if len(mu) else np.array([]);nv=ne[name,panel][1][take];by[name,panel]=(available,ev,nv)
                    if len(ev):results.append({'baseline':label,'scale':q,'panel':panel,'orientation':name,**summary_difference(nv,ev),'descriptive_only':True})
                    else:results.append({'baseline':label,'scale':q,'panel':panel,'orientation':name,'cells':0,'status':'OUTSIDE_MATCHED_SUPPORT','descriptive_only':True})
                    cells.extend({'baseline':label,'scale':q,'panel':panel,'orientation':name,'cell':int(n),'native_error':float(a),'baseline_error':float(b),'difference':float(a-b)} for n,a,b in zip(available,nv,ev))
            for ov,a,b in OVERLAPS:
                for panel in m:
                    aa,bb=by[a,panel],by[b,panel];nodes=np.intersect1d(np.intersect1d(aa[0],bb[0]),o[ov,q]);la=np.searchsorted(aa[0],nodes);lb=np.searchsorted(bb[0],nodes)
                    if len(nodes):results.append({'baseline':label,'scale':q,'panel':panel,'orientation':ov,**summary_difference((aa[2][la]+bb[2][lb])/2,(aa[1][la]+bb[1][lb])/2),'descriptive_only':True})
                    else:results.append({'baseline':label,'scale':q,'panel':panel,'orientation':ov,'cells':0,'status':'OUTSIDE_MATCHED_SUPPORT','descriptive_only':True})
    write('baselines/comparison.csv',results);write('baselines/support.csv',support);write('baselines/cellwise.csv.gz',cells);write('baselines/mnn_k5/native_support.csv',native_mnn)
    for label in tuple(embs)+('mnn_k5','mnn_fixed_size'):write(f'baselines/{label}/summary.csv',[r for r in results if r['baseline']==label])
    choices=[]
    for label in tuple(embs)+('mnn_k5','mnn_fixed_size'):
        for panel in m:
            rr=[r for r in results if r['baseline']==label and r['scale']=='q25' and r['panel']==panel and r['orientation'].startswith('O')];v=[r.get('paired_difference_mean') for r in rr]
            label0='OUTSIDE_CURRENT_AUTHORITY' if any(x is None for x in v) or len(v)!=3 else 'NATIVE_LOWER' if all(x<0 for x in v) else 'BASELINE_LOWER' if all(x>0 for x in v) else 'MIXED'
            choices.append({'baseline':label,'panel':panel,'overlap_differences':v,'classification':label0})
    js('baselines/decision.json',{'comparisons':choices,'population_p_values':False,'multiplicity_matched_mnn_completed':True,'biological_n':3});print('Baselines complete',choices,flush=True)
def representations():
    freeze();data={}
    for v in ('hash','hash_primary','pca32','pca64'):
        r,z,ix,c,g,o,m=arrays(v);ee,dd=errors(observed_means(g,'q25',m),g,o,'q25',m);data[v]=(c,g,o,ee,dd)
    objects=[];corr=[];directions=[]
    for ref in ('hash_primary','hash'):
        ca,ga,oa,ea,da=data[ref]
        for v in ('pca32','pca64'):
            cb,gb,ob,eb,db=data[v];meta={'reference':ref,'comparison':v,'like_for_like':ref=='hash_primary'}
            for p in PAIRS:
                for obj,i in [('source',0),('target',1)]:objects.append({**meta,'pair':p,'object':obj,**compare(ca[p][i],cb[p][i])})
                objects.append({**meta,'pair':p,'object':'edge',**compare(ca[p][0]*3704+ca[p][1],cb[p][0]*3704+cb[p][1])})
            for name in ga:
                objects.append({**meta,'orientation':name,'object':'domain',**compare(ga[name]['q25']['nodes'],gb[name]['q25']['nodes'])})
                for panel in ('state','heldout'):corr.append({**meta,'orientation':name,'panel':panel,'metric':'source_error',**correlation(ea[name,panel],eb[name,panel])})
            for name,_,_ in OVERLAPS:
                objects.append({**meta,'overlap':name,'object':'overlap',**compare(oa[name,'q25'],ob[name,'q25'])})
                for panel in ('state','heldout'):corr.append({**meta,'overlap':name,'panel':panel,'metric':'defect',**correlation(da[name,panel],db[name,panel])})
    for v in ('hash','pca32','pca64','target_pca32'):
        folder='matched_coherence/q25' if v=='hash' else 'representation/'+v
        rows0=read(ROOT/folder/'component_statistics.csv')
        for panel in ('state','heldout'):
            rr=[r for r in rows0 if r['panel']==panel and r['metric']=='defect'];directions.append({'representation':v,'panel':panel,'positive_component_signs':all(float(r['absolute_effect'])>0 for r in rr),
                'componentwise_supported':all(float(r['maxT_tail'])<=.05 and float(r['absolute_effect'])>0 for r in rr),'reductions':'|'.join(r['relative_effect'] for r in rr)})
    write('representation/object_identity.csv',objects);write('representation/cellwise_correlations.csv',corr);write('representation/effect_direction.csv',directions)
    for v in ('pca32','pca64'):write('representation/hash_vs_'+v+'.csv',[r for r in objects if r['comparison']==v])
    js('representation/decision.json',{'direction':directions,'all_objects_equal':all(r['Jaccard']==1 for r in objects),'object_identity_threshold_invented':False,'wording':'persists across representations when component signs remain positive, not independent replication'})
    print('Representation comparison complete',directions,flush=True)
def biology():
    freeze();full=rows(True);r,z,ix,c,g,o,m=arrays();f=np.load(P4/'geometry/hash/pair_cores.npz');oldc={p:(f[p+'_left'],f[p+'_right']) for p in PAIRS};oldz=z;og,oo=geometry(oldz,ix,oldc);records=[]
    for mode,cc,ov in [('v4',oldc,oo),('isolated',c,o)]:
        objects={}
        for p,(a,b) in ENDS.items():
            objects[p+'_source']=np.unique(cc[p][0]);objects[p+'_target']=np.unique(cc[p][1]);objects[p]=np.unique(np.r_[cc[p][0],cc[p][1]])
        for name,_,_ in OVERLAPS:objects[name]=ov[name,'q25']
        for name,nodes in objects.items():
            for field in ('cellType','germLayer','state','phase'):
                counts=Counter(full[int(i)][field] for i in nodes)
                for category in sorted({x[field] for x in full}):records.append({'normalization':mode,'object':name,'field':field,'category':category,'count':counts[category],'cells':len(nodes),'fraction':counts[category]/len(nodes) if len(nodes) else 0.,'descriptive_only':True})
    write('biology/landmark_annotation_summary.csv',[r for r in records if r['object'].startswith('P')]);write('biology/overlap_annotation_summary.csv',[r for r in records if r['object'].startswith('O')]);summary=[]
    for name in sorted({r['object'] for r in records}):
        for field in ('cellType','germLayer','state','phase'):
            whole=[r for r in records if r['normalization']=='v4' and r['object']==name and r['field']==field];isolated=[r for r in records if r['normalization']=='isolated' and r['object']==name and r['field']==field]
            a=max(whole,key=lambda r:(r['fraction'],r['category']));b=max(isolated,key=lambda r:(r['fraction'],r['category']));summary.append({'object':name,'field':field,'v4_dominant':a['category'],'v4_fraction':a['fraction'],'isolated_dominant':b['category'],'isolated_fraction':b['fraction'],'dominant_preserved':a['category']==b['category']})
    write('biology/object_composition_comparison.csv',summary);print('Descriptive biology complete',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['baselines','representations','biology']);a=p.parse_args();globals()[a.action]()

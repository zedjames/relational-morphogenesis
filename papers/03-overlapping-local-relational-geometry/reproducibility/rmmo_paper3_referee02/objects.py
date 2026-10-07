"""Competitive baselines, like-for-like object recurrence, descriptive biology."""
from collections import Counter
import itertools
import numpy as np
from scipy.stats import spearmanr
from engine import *
from ownership import arrays,observed_means

def errors(values,geo,overlaps,scale,matrices):
    e={};d={}
    for name,v in values.items():
        nodes=geo[name][scale]['nodes']
        for panel,m in matrices.items():
            mu=unit(v[panel]) if panel=='state' else v[panel]
            e[name,panel]=(nodes,np.linalg.norm(m[nodes]-mu,axis=1))
    for o,a,b in OVERLAPS:
        nodes=overlaps[o,scale]
        for panel in matrices:
            va=values[a][panel][np.searchsorted(geo[a][scale]['nodes'],nodes)]
            vb=values[b][panel][np.searchsorted(geo[b][scale]['nodes'],nodes)]
            if panel=='state':va,vb=unit(va),unit(vb)
            d[o,panel]=(nodes,np.linalg.norm(va-vb,axis=1))
    return e,d

def summary_difference(obs,base):
    diff=obs-base
    return {'cells':len(diff),'observed_mean':float(obs.mean()),'baseline_mean':float(base.mean()),
        'paired_difference_mean':float(diff.mean()),'paired_difference_median':float(np.median(diff)),
        'observed_median':float(np.median(obs)),'baseline_median':float(np.median(base))}

def baseline():
    freeze();rows,z,ix,core,geo,overlaps,matrices=arrays();stats=[];defects=[];support=[];cells=[]
    embeddings={'matched_knn':z}
    for k in (32,64):
        with np.load(PARENT/f'representation_families/pca{k}/representation.npz') as f:embeddings[f'pca{k}_knn']=f['embedding'].copy()
    distances={name:pair_distances(emb,ix)[0] for name,emb in embeddings.items()}
    for scale in SCALES:
        observed=observed_means(geo,scale,matrices);err,de=errors(observed,geo,overlaps,scale,matrices)
        for label,ds in distances.items():
            selected={};values={}
            for name,e in geo.items():
                p='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))))
                distance=ds[p] if e['source']<e['target'] else ds[p].T
                tar,mem=fiber_matrix(e,scale);sizes=mem.sum(1)
                order=np.argsort(distance[np.searchsorted(ix[e['source']],e[scale]['nodes'])],axis=1,kind='stable')
                selected[name]=[ix[e['target']][order[i,:k]] for i,k in enumerate(sizes)]
                values[name]={panel:np.asarray([m[t].mean(0) for t in selected[name]]) for panel,m in matrices.items()}
            be,bd=errors(values,geo,overlaps,scale,matrices)
            for name,e in geo.items():
                support.append({'baseline':label,'scale':scale,'orientation':name,'native_support':len(e[scale]['nodes']),
                    'common_support':len(e[scale]['nodes']),'mean_multiplicity':float(np.mean([len(t) for t in selected[name]])),
                    'source_support_and_size_matched':True})
                for panel in matrices:
                    nodes,obs=err[name,panel];bv=be[name,panel][1]
                    stats.append({'baseline':label,'scale':scale,'panel':panel,'orientation':name,**summary_difference(obs,bv)})
                    if label=='matched_knn':
                        cells.extend({'scale':scale,'orientation':name,'panel':panel,'source_node':int(n),'observed_error':float(a),'baseline_error':float(b)} for n,a,b in zip(nodes,obs,bv))
            for o,a,b in OVERLAPS:
                nodes=overlaps[o,scale]
                for panel in matrices:
                    defects.append({'baseline':label,'scale':scale,'panel':panel,'overlap':o,**summary_difference(de[o,panel][1],bd[o,panel][1])})
                    ae=(err[a,panel][1][np.searchsorted(err[a,panel][0],nodes)]+err[b,panel][1][np.searchsorted(err[b,panel][0],nodes)])/2
                    ab=(be[a,panel][1][np.searchsorted(be[a,panel][0],nodes)]+be[b,panel][1][np.searchsorted(be[b,panel][0],nodes)])/2
                    stats.append({'baseline':label,'scale':scale,'panel':panel,'orientation':o,**summary_difference(ae,ab)})
        # MNN5 uses native support, with explicitly common-support comparisons.
        ds=distances['matched_knn'];native={};mnnvalues={};mnnsizes={}
        for name,e in geo.items():
            p='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))))
            distance=ds[p] if e['source']<e['target'] else ds[p].T
            f=np.argsort(distance,axis=1,kind='stable')[:,:5]
            rev=np.argsort(distance,axis=0,kind='stable')[:5,:]
            member=np.zeros_like(distance,dtype=bool);member[np.arange(len(distance))[:,None],f]=True
            reverse=np.zeros_like(distance,dtype=bool);reverse[rev,np.arange(distance.shape[1])[None,:]]=True
            member &= reverse
            available=member.any(1);native[name]=ix[e['source']][available]
            sizes=member[available].sum(1);mnnsizes[name]=sizes
            mnnvalues[name]={panel:(member[available].astype(float)@m[ix[e['target']]])/sizes[:,None] for panel,m in matrices.items()}
            common=np.intersect1d(e[scale]['nodes'],native[name]);sel=np.searchsorted(native[name],common)
            support.append({'baseline':'mnn_k5','scale':scale,'orientation':name,'native_support':len(native[name]),
                'common_support':len(common),'mean_multiplicity':float(sizes.mean()),'source_support_and_size_matched':False})
            for panel,m in matrices.items():
                mm=mnnvalues[name][panel];mm=unit(mm) if panel=='state' else mm
                bv=np.linalg.norm(m[common]-mm[sel],axis=1);obs=err[name,panel][1][np.searchsorted(err[name,panel][0],common)]
                stats.append({'baseline':'mnn_k5','scale':scale,'panel':panel,'orientation':name,**summary_difference(obs,bv)})
        for o,a,b in OVERLAPS:
            common=np.intersect1d(np.intersect1d(overlaps[o,scale],native[a]),native[b])
            for panel,m in matrices.items():
                va=mnnvalues[a][panel][np.searchsorted(native[a],common)];vb=mnnvalues[b][panel][np.searchsorted(native[b],common)]
                if panel=='state':va,vb=unit(va),unit(vb)
                bv=np.linalg.norm(va-vb,axis=1);obs=de[o,panel][1][np.searchsorted(de[o,panel][0],common)]
                defects.append({'baseline':'mnn_k5','scale':scale,'panel':panel,'overlap':o,**summary_difference(obs,bv)})
                oa=err[a,panel][1][np.searchsorted(err[a,panel][0],common)];ob=err[b,panel][1][np.searchsorted(err[b,panel][0],common)]
                eb=(np.linalg.norm(m[common]-va,axis=1)+np.linalg.norm(m[common]-vb,axis=1))/2
                stats.append({'baseline':'mnn_k5','scale':scale,'panel':panel,'orientation':o,**summary_difference((oa+ob)/2,eb)})
    write('competitive_baselines/source_fidelity.csv',stats);write('competitive_baselines/heldout_fidelity.csv',[r for r in stats if r['panel']=='heldout'])
    write('competitive_baselines/overlap_defect.csv',defects);write('competitive_baselines/support_summary.csv',support)
    write('competitive_baselines/pairwise_comparison.csv',stats);write('competitive_baselines/matched_knn/cellwise_errors.csv.gz',cells)
    decisions=[]
    for label in ('matched_knn','mnn_k5','pca32_knn','pca64_knn'):
        for panel in matrices:
            vv=[r for r in stats if r['baseline']==label and r['scale']=='q25' and r['panel']==panel and r['orientation'].startswith('O')]
            dif=[r['paired_difference_mean'] for r in vv]
            decision='LANDMARK_FIDELITY_ADVANTAGE' if all(x < -1e-8 for x in dif) else 'SIMPLE_BASELINE_OUTPERFORMS' if all(x>1e-8 for x in dif) else 'EQUIVALENT_TO_SIMPLE_LOCAL_GEOMETRY' if all(abs(x)<=1e-8 for x in dif) else 'MIXED_BY_EMBRYO_OR_REPRESENTATION'
            decisions.append({'baseline':label,'panel':panel,'paired_differences':dif,'decision':decision})
        folder='pca_knn' if label.startswith('pca') else label
        write(f'competitive_baselines/{folder}/{label}_summary.csv',[r for r in stats if r['baseline']==label])
    js('competitive_baselines/decision.json',{'comparisons':decisions,'biological_n':3,'cellwise_population_p_values':False})
    print('Competitive baselines',decisions,flush=True)

def jaccard(a,b):
    a,b=set(a),set(b);return len(a&b)/len(a|b) if a|b else 1.
def correlation(a,b):
    nodes=np.intersect1d(a[0],b[0]);aa=a[1][np.searchsorted(a[0],nodes)];bb=b[1][np.searchsorted(b[0],nodes)]
    rho=spearmanr(aa,bb).statistic if len(nodes)>2 and np.std(aa)>0 and np.std(bb)>0 else None
    return {'common_cells':len(nodes),'spearman_rho':float(rho) if rho is not None and np.isfinite(rho) else None}

def representation():
    freeze();rows,z,ix,stable,_,_,_=arrays();bank=npz_parent(15,'pair_native/lane_relations.npz')
    family={'hash_primary_q50':(z,{p:(bank['primary_'+p+'_left'],bank['primary_'+p+'_right']) for p in PAIRS}),
            'hash_stable':(z,stable)}
    for k in (32,64):
        with np.load(PARENT/f'representation_families/pca{k}/representation.npz') as f:emb=f['embedding'].copy()
        cc={}
        for p,(a,b) in ENDS.items():
            d=chord(emb[ix[a]],emb[ix[b]]);eps=np.quantile(np.r_[d.min(1),d.min(0)],.5);left,right=np.where(d<=eps)
            cc[p]=(ix[a][left],ix[b][right])
        family[f'pca{k}']=(emb,cc)
    objects={};endpoint=[];domain=[];overlap=[];sourcecorr=[];defectcorr=[]
    for label,(emb,cc) in family.items():
        geo,oo=geometry(emb,ix,cc);matrices={'state':emb,'heldout':load_target()}
        vv=observed_means(geo,'q25',matrices);ee,dd=errors(vv,geo,oo,'q25',matrices)
        objects[label]=(geo,oo,ee,dd)
    for reference in ('hash_primary_q50','hash_stable'):
        for label in ('pca32','pca64'):
            meta={'reference':reference,'comparison':label,'like_for_like':reference=='hash_primary_q50'}
            for p in PAIRS:
                aa,bb=family[reference][1][p],family[label][1][p]
                endpoint.append({**meta,'pair':p,'source_endpoint_jaccard':jaccard(aa[0],bb[0]),
                    'target_endpoint_jaccard':jaccard(aa[1],bb[1]),'edge_jaccard':jaccard(zip(*aa),zip(*bb))})
            ga,oa,ea,da=objects[reference];gb,ob,eb,db=objects[label]
            for name in ga:
                domain.append({**meta,'orientation':name,'domain_jaccard':jaccard(ga[name]['q25']['nodes'],gb[name]['q25']['nodes'])})
                for panel in ('state','heldout'):sourcecorr.append({**meta,'orientation':name,'panel':panel,**correlation(ea[name,panel],eb[name,panel])})
            for o,_,_ in OVERLAPS:
                overlap.append({**meta,'overlap':o,'overlap_cell_jaccard':jaccard(oa[o,'q25'],ob[o,'q25'])})
                for panel in ('state','heldout'):defectcorr.append({**meta,'overlap':o,'panel':panel,**correlation(da[o,panel],db[o,panel])})
    write('cross_representation/landmark_overlap.csv',endpoint);write('cross_representation/domain_overlap.csv',domain)
    write('cross_representation/overlap_cell_overlap.csv',overlap);write('cross_representation/source_error_correlations.csv',sourcecorr)
    write('cross_representation/overlap_defect_correlations.csv',defectcorr)
    summary=[]
    matched=read(ROOT/'family_statistics/component_statistics.csv')
    for label in ('pca32','pca64'):
        rr=[r for r in endpoint if r['comparison']==label and r['like_for_like']]
        directional=[r for r in matched if r['authority']==label and r['bank']=='matched' and r['scale']=='q25' and r['metric']=='defect']
        assert len(directional)==6
        summary.append({'comparison':label,'endpoint_recurrence':all(r[k]>=.5 for r in rr for k in ('source_endpoint_jaccard','target_endpoint_jaccard')),
            'edge_recurrence':all(r['edge_jaccard']>=.5 for r in rr),
            'domain_recurrence':all(r['domain_jaccard']>=.5 for r in domain if r['comparison']==label and r['like_for_like']),
            'overlap_cell_recurrence':all(r['overlap_cell_jaccard']>=.5 for r in overlap if r['comparison']==label and r['like_for_like']),
            'source_error_cellwise_statistic_recurrence':all(r['spearman_rho'] is not None and r['spearman_rho']>=.5 for r in sourcecorr if r['comparison']==label and r['like_for_like']),
            'overlap_defect_cellwise_statistic_recurrence':all(r['spearman_rho'] is not None and r['spearman_rho']>=.5 for r in defectcorr if r['comparison']==label and r['like_for_like']),
            'matched_coherence_effect_direction_recurrence':all(float(r['absolute_effect'])>0 for r in directional),
            'matched_coherence_10_percent_gate_recurrence':all(float(r['relative_effect'])>=.10 for r in directional if r['panel']=='realization')})
    write('cross_representation/object_recurrence_summary.csv',summary)
    js('cross_representation/decision.json',{'objects':summary,'direction_authority':'immutable Referee01 cross-family direction, not full effect-size-gate preservation','threshold':.5,'threshold_is_descriptive':True})
    print('Representation object recurrence',summary,flush=True)

def biology():
    freeze();rows,z,ix,core,geo,overlaps,matrices=arrays();ends=[];counts=[];enrichment=[];oc=[];oe=[]
    fields=('cellType','germLayer','state','phase','sourceLineage')
    def describe(nodes,rep,metadata,countout,enrichout):
        for field in fields:
            background=Counter(rows[int(n)][field] for n in ix[rep]);actual=Counter(rows[int(n)][field] for n in nodes)
            for category,bn in sorted(background.items()):
                n=actual[category];f=n/len(nodes) if len(nodes) else 0.;bf=bn/len(ix[rep])
                r={**metadata,'field':field,'category':category,'count':n,'fraction':f,'background_count':bn,
                    'background_fraction':bf,'observed_background_ratio':f/bf,'proportion_difference':f-bf}
                countout.append(r);enrichout.append(r)
    for p,(a,b) in ENDS.items():
        for side,rep,nodes in [('left',a,np.unique(core[p][0])),('right',b,np.unique(core[p][1]))]:
            for n in nodes:ends.append({'source_node':int(n),'replicate':old.REPS[rep],'pair':p,'side':side,**{f:rows[int(n)][f] for f in fields}})
            describe(nodes,rep,{'pair':p,'side':side,'replicate':old.REPS[rep]},counts,enrichment)
    for o,a,b in OVERLAPS:describe(overlaps[o,'q25'],geo[a]['source'],{'overlap':o,'replicate':old.REPS[geo[a]['source']]},oc,oe)
    write('biology/landmark_endpoints.csv',ends);write('biology/landmark_annotation_counts.csv',counts)
    write('biology/landmark_annotation_enrichment.csv',enrichment);write('biology/overlap_annotation_counts.csv',oc)
    write('biology/overlap_annotation_enrichment.csv',oe)
    cells=read(ROOT/'competitive_baselines/matched_knn/cellwise_errors.csv.gz')
    vv=observed_means(geo,'q25',matrices);err,_=errors(vv,geo,overlaps,'q25',matrices)
    # Independent exact stored degree-only draws, not an invented matched control.
    with np.load(PARENT/'representation_families/hash_family_reference/degree_only_uniform_stub_rejection_rewired_cores.npz') as f:draws={p:f[p+'_target'].copy() for p in PAIRS}
    breaksum={(name,p):np.zeros(len(geo[name]['q25']['nodes'])) for name in geo for p in matrices}
    for b in range(B):
        new={p:(core[p][0],draws[p][b]) for p in PAIRS}
        ev,_=errors(observed_means(geo,'q25',matrices,new),geo,overlaps,'q25',matrices)
        for key in breaksum:breaksum[key]+=ev[key][1]/B
    for field,filename in [('germLayer','germ_layer'),('state','S0'),('cellType','cell_type')]:
        out=[]
        for name in geo:
            for panel in matrices:
                rr=[r for r in cells if r['scale']=='q25' and r['orientation']==name and r['panel']==panel]
                groups={rows[int(r['source_node'])][field] for r in rr}
                for category in sorted(groups):
                    ss=[r for r in rr if rows[int(r['source_node'])][field]==category]
                    if len(ss)<20:continue
                    nodes=np.asarray([int(r['source_node']) for r in ss]);loc=np.searchsorted(geo[name]['q25']['nodes'],nodes)
                    observed=float(np.mean([float(r['observed_error']) for r in ss]));base=float(np.mean([float(r['baseline_error']) for r in ss]));broken=float(breaksum[name,panel][loc].mean())
                    out.append({'orientation':name,'panel':panel,'category':category,'source_count':len(ss),'observed_mean_error':observed,
                        'matched_knn_mean_error':base,'pair_native_minus_baseline':observed-base,
                        'degree_break_mean_error':broken,'correspondence_break_effect':broken-observed,
                        'correspondence_break_reference':'degree_only_uniform_stub_rejection','descriptive_only':True})
        write(f'biology/source_fidelity_by_{filename}.csv',out)
    text_file('biology/interpretation.md','# Biological description\n\nAll annotation categories, including zero-count categories, are supplied against each embryo background. Labels are expression-derived and are never independent validation. Compartment fidelity summaries require at least20 cells and describe this measured carrier only. Their degree-break effects use the immutable1024 exact degree-only draws. Neither P13 composition nor any enrichment ratio establishes a causal explanation. Biological n=3; no population generalization.')
    print('Descriptive biological tables complete',flush=True)

if __name__=='__main__':
    import sys
    {'baseline':baseline,'representation':representation,'biology':biology}[sys.argv[1]]()

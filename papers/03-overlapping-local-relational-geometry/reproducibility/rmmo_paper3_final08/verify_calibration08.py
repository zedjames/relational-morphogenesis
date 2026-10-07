"""Independent scalar ranking and sparse-incidence calibration replay; never writes banks."""
import argparse, statistics
import scipy.sparse as sp
from io08 import *
from publication_stats08 import families, calibration_source
from calibration08 import inputs
from banks08 import Sampler, context, check, seed, rng, p

def scalar_ranks():
    constants=read(ROOT/'statistical_authority/component_calibration_final.csv')
    rows=read(ROOT/'statistical_authority/final_component_statistics.csv')
    comparisons=read(ROOT/'statistical_authority/final_maxT_comparison.csv')
    family=read(ROOT/'statistical_authority/final_family_statistics.csv')
    events=read(ROOT/'statistical_authority/final_family_decision_comparison.csv')
    select=lambda rr,f:[r for r in rr if all(r[k]==str(f[k]) for k in ('bank','panel','metric'))]
    for f in families():
        _,old,v,obs=inputs(f); path,offset,_=calibration_source(f)
        cv=np.load(path)['statistics'][:,offset:offset+3];cc=select(constants,f);rr=select(rows,f);cmp=select(comparisons,f)
        mu=[float(r['mu_cal']) for r in cc];sd=[float(r['sigma_cal']) for r in cc];s=f['direction']
        for i in range(3):
            np.testing.assert_allclose(mu[i],statistics.mean(map(float,cv[:,i])),rtol=1e-12,atol=1e-12)
            np.testing.assert_allclose(sd[i],statistics.stdev(map(float,cv[:,i])),rtol=1e-12,atol=1e-12)
        zo=[s*(float(obs[i])-mu[i])/(sd[i]+1e-8) for i in range(3)]
        zz=[[s*(float(x[i])-mu[i])/(sd[i]+1e-8) for i in range(3)] for x in v]
        maxima=[max(x) for x in zz]
        for i in range(3):
            tail=(1+sum(x>=zo[i] for x in maxima))/1025
            raw=(1+sum(s*float(x[i])>=s*float(obs[i]) for x in v))/1025
            assert tail==float(rr[i]['maxT_tail'])==float(cmp[i]['calibrated_maxT_tail'])
            assert raw==float(rr[i]['conditional_tail'])==float(old[i]['conditional_tail'])
            assert float(cmp[i]['original_maxT_tail'])==float(old[i]['maxT_tail'])
            assert cmp[i]['decision_changed']==str((tail<=.05)!=(float(old[i]['maxT_tail'])<=.05))
            np.testing.assert_allclose(float(rr[i]['Z']),zo[i],rtol=1e-12,atol=1e-12)
        ff=select(family,f)[0]
        for name,reduce in [('sum',sum),('min',min)]:
            tail=(1+sum(reduce(x)>=reduce(zo) for x in zz))/1025
            assert tail==float(ff[name+'_tail'])
            event=next(x for x in select(events,f) if x['statistic']==name)
            assert tail==float(event['calibrated_tail'])
            assert event['decision_changed']==str((tail<=.05)!=(float(event['original_tail'])<=.05))
    assert len(rows)==108 and len(family)==36 and len(events)==72
    receipt=json.loads((ROOT/'freeze/component_constants_freeze.json').read_text())
    assert sha(ROOT/'statistical_authority/component_calibration_final.csv')==receipt['constants_sha256']
    js('completion/final_scalar_validation.json',dict(passed=True,families=36,components=108,family_events=72,raw_tails_unchanged=True,independent_scalar_rank_reproduction=True,sample_sd_ddof=1,eta=1e-8))
    print('PASS independently reproduced all 108 component ranks and 72 family ranks',flush=True)

def independent_scores(g,o,m,values,q,metric):
    result=[]
    for panel in ('state','heldout'):
        v={name:p.unit(x[panel]) if panel=='state' else x[panel] for name,x in values.items()}
        for name,a,b in OVERLAPS:
            nodes=o[name,q];first=v[a][np.searchsorted(g[a][q]['nodes'],nodes)];second=v[b][np.searchsorted(g[b][q]['nodes'],nodes)]
            if metric=='defect': distances=np.sqrt(np.sum((first-second)**2,axis=1))
            else: distances=(np.sqrt(np.sum((m[panel][nodes]-first)**2,axis=1))+np.sqrt(np.sum((m[panel][nodes]-second)**2,axis=1)))/2
            result.append(float(distances.mean()))
    return np.asarray(result)

def matched_replay():
    plan=json.loads((ROOT/'freeze/missing_calibration_plan.json').read_text());count=0;seeds=set()
    for spec in plan['banks']:
        if spec['metric']!='defect':continue
        bank,q=spec['bank'],spec['scale'];r,z,ix,c,g,o,m=p.arrays(spec['variant']);ss={n:Sampler(2,r,z,ix,e,q) for n,e in g.items()};contexts={n:context(s) for n,s in ss.items()};path=ROOT/'independent_calibration_bank'/bank
        log={(int(x['draw']),x['orientation']):x for x in read(path/'graph_audit.csv')};stored=np.load(path/'bank.npz');assert len(log)==6144
        for b in range(B):
            values={}
            for name,s in ss.items():
                label=f'matched|{bank}|{b}|{name}';s.draw(rng(label),{});check(s,contexts[name]);codes=s.last_graph_codes;record=log[b,name]
                assert s.last_graph_hash==record['graph_hash'] and seed(label)==int(record['seed'])
                assert record['replacement']=='False' and record['fallback']=='0' and record['NOT_INFERENTIAL']=='True';seeds.add(seed(label))
                incidence=sp.csr_matrix((np.ones(len(codes)),(np.searchsorted(s.nodes,codes//3704),codes%3704)),shape=(len(s.nodes),3704))
                values[name]={panel:np.asarray(incidence@matrix)/s.size[:,None] for panel,matrix in m.items()};count+=1
            np.testing.assert_allclose(independent_scores(g,o,m,values,q,'defect'),stored['statistics'][b],rtol=1e-12,atol=1e-12)
            if b%256==0:print('Independent matched incidence replay',bank,b,'of',B,flush=True)
        observed={name:{panel:np.asarray(sp.csr_matrix(s.observed.astype(float))@matrix[s.targets])/s.size[:,None] for panel,matrix in m.items()} for name,s in ss.items()}
        np.testing.assert_allclose(independent_scores(g,o,m,observed,q,'defect'),stored['observed'],rtol=1e-12,atol=1e-12)
    assert count==30720 and len(seeds)==30720
    js('completion/matched_calibration_replay.json',dict(passed=True,banks=5,draws_per_bank=1024,exact_graphs=count,distinct_seeds=len(seeds),independent_sparse_incidence_scores=True,all_observed_scores_verified=True))
    print('PASS five matched calibration banks / 30720 graphs',flush=True)

def degree_replay():
    r,z,ix,c,g,o,m=p.arrays();base=ROOT/'independent_calibration_bank/source_fidelity';graphs=np.load(base/'degree_only_graphs.npz');logs={(int(x['draw']),x['pair']):x for x in read(base/'graph_audit.csv')};banks={q:np.load(base/q/'bank.npz') for q in ('q10','q50')};count=0
    for b in range(B):
        new={}
        for name,(left,right) in c.items():
            label=f'degree_missing|{name}|{b}';random=rng(label);trials=0
            while True:
                t=random.permutation(right);trials+=1
                if len(np.unique(left*3704+t))==len(left):break
            record=logs[b,name];np.testing.assert_array_equal(t,graphs[name][b]);np.testing.assert_array_equal(np.sort(t),np.sort(right))
            assert seed(label)==int(record['seed']) and trials==int(record['stub_trials']) and p.digest_graph(left,t)==record['graph_hash'];new[name]=(left,t);count+=1
        for q in banks:
            values={}
            for name,e in g.items():
                pair='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))));left,right=new[pair];se,te=(left,right) if e['source']<e['target'] else (right,left)
                adjacency=sp.csr_matrix((np.ones(len(se)),(np.searchsorted(e['anchors'],se),te)),shape=(len(e['anchors']),3704))
                incidence=(sp.csr_matrix(e[q]['active'].astype(float))@adjacency).astype(bool).astype(float);sizes=np.asarray(incidence.sum(1)).ravel();assert (sizes>0).all()
                values[name]={panel:np.asarray(incidence@matrix)/sizes[:,None] for panel,matrix in m.items()}
            np.testing.assert_allclose(independent_scores(g,o,m,values,q,'source_error'),banks[q]['statistics'][b],rtol=1e-12,atol=1e-12)
        if b%256==0:print('Independent degree incidence replay',b,'of',B,flush=True)
    assert count==3072
    js('completion/degree_calibration_replay.json',dict(passed=True,banks=2,shared_pair_graphs=count,distinct_seeds=count,draws_per_bank=1024,exact_degree_margins_simplicity_and_endpoints=True,independent_sparse_incidence_scores=True))
    print('PASS two degree calibration banks / 3072 shared graphs',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['scalar_ranks','matched_replay','degree_replay']);globals()[parser.parse_args().action]()

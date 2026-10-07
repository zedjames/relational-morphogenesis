"""Direct IID retained laws; no switch chains or rescues of closed claims."""
import argparse
from common05 import *
def edges():
    freeze();r,z,ix,c,g,o,m=arrays();_,_,bins=pair_distances(z,ix);states=sorted({x['state'] for x in r});sc={s:i for i,s in enumerate(states)};st=np.asarray([sc[x['state']] for x in r]);ns=len(states)
    ctx={};obs=[];audit=[]
    for p,(a,b) in ENDS.items():
        left,right=c[p];j=(st[ix[a]][:,None]*ns+st[ix[b]][None,:])*10+bins[p]
        k=(st[left]*ns+st[right])*10+bins[p][np.searchsorted(ix[a],left),np.searchsorted(ix[b],right)];fields={}
        for key,n in zip(*np.unique(k,return_counts=True)):
            l,t=np.where(j==key);fields[int(key)]=(ix[a][l],ix[b][t],int(n));assert n<=len(l)
            audit.append({'pair':p,'stratum':int(key),'observed_count':int(n),'pool_edges':len(l),'state_only_matching':True})
        ctx[p]=fields;obs.append(float(np.linalg.norm(m['heldout'][left]-m['heldout'][right],axis=1).mean()))
    null=[];v=[];graphs=[]
    for d in range(B):
        vv=[]
        for p,fields in ctx.items():
            random=rng(f'edge|{p}|{d}');ll=[];tt=[]
            for k,(l,t,n) in fields.items():
                pick=random.choice(len(l),size=n,replace=False);ll.extend(l[pick]);tt.extend(t[pick])
            ll,tt=np.asarray(ll),np.asarray(tt);assert len(set(zip(ll,tt)))==len(ll)
            a,b=ENDS[p];actual=(st[ll]*ns+st[tt])*10+bins[p][np.searchsorted(ix[a],ll),np.searchsorted(ix[b],tt)]
            assert Counter(actual)==Counter({k:v[2] for k,v in fields.items()})
            val=float(np.linalg.norm(m['heldout'][ll]-m['heldout'][tt],axis=1).mean());vv.append(val)
            null.append({'draw':d,'pair':p,'mean_distance':val});graphs.append({'draw':d,'pair':p,'seed':seed(f'edge|{p}|{d}'),'graph_hash':digest_graph(ll,tt),'constraints_validated':True})
        v.append(vv)
    save('edge_agreement/bank.npz',observed=np.asarray(obs),statistics=np.asarray(v));write('edge_agreement/nulls_1024.csv',null);write('edge_agreement/graph_audit.csv',graphs);write('edge_agreement/strata.csv',audit)
    write('edge_agreement/observed.csv',[{'pair':p,'edges':len(c[p][0]),'mean_distance':x} for p,x in zip(PAIRS,obs)]);d=summarize('edge_agreement',obs,v,PAIRS)
    print('Edge agreement IID complete',[(p,x['relative_effect'],x['maxT_tail']) for p,x in zip(PAIRS,d['components'])],flush=True)
def spatial():
    freeze();r,z,ix,c,g,o,m=arrays();ds={};den={};nn={};contexts={};spec=[];exports={}
    for rep,nodes in enumerate(ix):
        d=chord(z[nodes],z[nodes]);ds[rep]=d;v=d.copy();np.fill_diagonal(v,np.inf)
        ten=np.partition(v,9,axis=1)[:,9];one=v.min(1)
        den[rep]=np.searchsorted(np.quantile(ten,[.25,.5,.75]),ten,side='right');nn[rep]=np.searchsorted(np.quantile(one,np.arange(.1,1,.1)),one,side='right')
    for name,e in g.items():
        rep=e['source'];nodes=ix[rep];obs=np.searchsorted(nodes,e['anchors']);center=z[e['anchors']].mean(0);rad=np.linalg.norm(z[nodes]-center,axis=1)
        rb=np.searchsorted(np.quantile(rad,np.arange(.1,1,.1)),rad,side='right')
        for ref in ('S1',):
            keys=[(r[int(v)]['state'],int(den[rep][i]))+((int(rb[i]),int(nn[rep][i])) if ref=='S1' else ()) for i,v in enumerate(nodes)]
            pools=defaultdict(list)
            for i,k in enumerate(keys):pools[k].append(i)
            hist=Counter(keys[int(i)] for i in obs);contexts[ref,name]=(keys,pools,hist)
            spec.extend({'reference':ref,'orientation':name,'stratum':str(k),'count':n,'pool_size':len(pools[k]),'fixed_precomputed':True} for k,n in hist.items())
            exports[ref+'_'+name+'_keys']=np.asarray([str(k) for k in keys]);exports[ref+'_'+name+'_nodes']=nodes
    save('spatial_overlap/fixed_cell_strata.npz',**exports);write('spatial_overlap/fixed_strata.csv',spec)
    js('spatial_overlap/S1_construction_audit.json',{'fixed_before_sampling':True,'strata_sha256':sha(ROOT/'spatial_overlap/fixed_cell_strata.npz'),
        'exact_wording':'Each S1 reference landmark system is formed by independent uniform sampling without replacement within fixed precomputed cell-level strata.','recorded_utc':utc(),'S2_rerun':False})
    for ref in ('S1',):
        vals=[];records=[];draws={name:[] for name in g};hashes=[]
        for b in range(B):
            domains={}
            for name,e in g.items():
                keys,pools,hist=contexts[ref,name];random=rng(f'spatial|{ref}|{b}|{name}')
                selected=np.concatenate([random.choice(pools[k],size=n,replace=False) for k,n in hist.items()]);assert len(np.unique(selected))==len(selected) and Counter(keys[int(i)] for i in selected)==hist
                nodes=ix[e['source']][selected];draws[name].append(nodes);rho=ds[e['source']][:,selected].min(1);rad=float(np.quantile(rho,.25));domains[name]=ix[e['source']][rho<=rad]
                hashes.append({'reference':ref,'draw':b,'orientation':name,'seed':seed(f'spatial|{ref}|{b}|{name}'),'landmark_hash':hashlib.sha256(np.sort(nodes).astype(np.int64).tobytes()).hexdigest(),'matching_validated':True,'radius':rad})
            vv=[len(np.intersect1d(domains[a],domains[b])) for _,a,b in OVERLAPS];vals.append(vv)
            records.extend({'reference':ref,'draw':b,'overlap':name,'size':n} for (name,_,_),n in zip(OVERLAPS,vv))
        folder='spatial_overlap/'+('original_reference' if ref=='original' else 'S1');observed=[len(o[name,'q25']) for name,_,_ in OVERLAPS]
        save(folder+'/landmark_draws.npz',**{name:np.asarray(v) for name,v in draws.items()});save(folder+'/bank.npz',observed=np.asarray(observed),statistics=np.asarray(vals));write(folder+'/nulls_1024.csv',records);write(folder+'/graph_audit.csv',hashes)
        d=summarize(folder,observed,vals,[name for name,_,_ in OVERLAPS],False,{'reference':ref});print('Spatial IID complete',ref,d['sum_tail'],flush=True)
    for filename in ('component_statistics.csv','family_statistics.csv'):write('spatial_overlap/'+filename,[x for rel in ('S1',) for x in read(ROOT/f'spatial_overlap/{rel}/{filename}')])
    js('spatial_overlap/decision.json',{ref:json.loads((ROOT/f'spatial_overlap/{rel}/decision.json').read_text()) for ref,rel in [('S1','S1')]})
def context(s):
    keys={k:i for i,k in enumerate(sorted(set(s.key.values())))};codes=s.nodes[np.where(s.observed)[0]]*3704+s.targets[np.where(s.observed)[1]]
    lookup=np.asarray([keys[s.key[int(t)]] if int(t) in s.key else -1 for t in range(3704)])
    expected=np.sort((codes//3704)*len(keys)+lookup[codes%3704]);return len(keys),lookup,expected
def check_incidence(s,ctx):
    nc,lookup,expected=ctx;codes=s.last_graph_codes;actual=np.sort((codes//3704)*nc+lookup[codes%3704]);assert np.array_equal(actual,expected)
def matched(variant='hash'):
    freeze();sv='hash' if variant=='target_pca32' else variant;r,z,ix,c,g,o,m=arrays(sv,variant);scales=SCALES if variant=='hash' else ['q25']
    for q in scales:
        samplers={n:Sampler(2,r,z,ix,e,q) for n,e in g.items()};contexts={n:context(s) for n,s in samplers.items()};bank=[];hashes=[]
        for b in range(B):
            values={}
            for name,s in samplers.items():
                values[name],used=s.draw(rng(f'matched|{variant}|{q}|{b}|{name}'),m);check_incidence(s,contexts[name])
                hashes.append({'draw':b,'orientation':name,'seed':seed(f'matched|{variant}|{q}|{b}|{name}'),'graph_hash':s.last_graph_hash,'matching_validated':True})
            bank.append(score(g,o,m,q,values=values))
            if b%256==0:print('Matched IID',variant,q,b,flush=True)
        rel='matched_coherence/'+q if variant=='hash' else 'representation/'+variant;obs=score(g,o,m,q);save(rel+'/bank.npz',observed=obs,statistics=np.asarray(bank));write(rel+'/graph_audit.csv',hashes);vector_summaries(rel,obs,np.asarray(bank))
    if variant=='hash':
        for name in ('component_statistics.csv','family_statistics.csv'):write('matched_coherence/'+name,[{'scale':q,**r} for q in SCALES for r in read(ROOT/f'matched_coherence/{q}/{name}')])
        data=read(ROOT/'matched_coherence/component_statistics.csv')
        for panel in ('state','heldout'):write(f'matched_coherence/{panel}_summary.csv',[r for r in data if r['panel']==panel]);
        js('matched_coherence/decision.json',{q:json.loads((ROOT/f'matched_coherence/{q}/decision.json').read_text()) for q in SCALES})
def degree():
    freeze();r,z,ix,c,g,o,m=arrays();banks={q:[] for q in SCALES};edges={p:[] for p in PAIRS};audit=[]
    for d in range(B):
        new={}
        for p,(left,right) in c.items():
            random=rng(f'degree|{p}|{d}');tries=0
            while True:
                shuffled=random.permutation(right);tries+=1
                if len(np.unique(left*3704+shuffled))==len(left):break
            assert np.array_equal(np.sort(right),np.sort(shuffled));assert set(shuffled)==set(right);new[p]=(left,shuffled);edges[p].append(shuffled)
            audit.append({'draw':d,'pair':p,'seed':seed(f'degree|{p}|{d}'),'stub_trials':tries,'graph_hash':digest_graph(left,shuffled),'simplicity_degrees_endpoints_validated':True})
        for q in SCALES:banks[q].append(score(g,o,m,q,edges=new))
        if d%256==0:print('Degree IID',d,flush=True)
    save('source_fidelity/degree_only_graphs.npz',**{p:np.asarray(v) for p,v in edges.items()});write('source_fidelity/graph_audit.csv',audit)
    for q,v in banks.items():
        obs=score(g,o,m,q);save(f'source_fidelity/{q}/bank.npz',observed=obs,statistics=np.asarray(v));vector_summaries('source_fidelity/'+q,obs,np.asarray(v))
    for name in ('component_statistics.csv','family_statistics.csv'):write('source_fidelity/'+name,[{'scale':q,**r} for q in SCALES for r in read(ROOT/f'source_fidelity/{q}/{name}')])
    write('source_fidelity/degree_only_nulls_1024.csv',[{'scale':q,'draw':b,'statistics':'|'.join(map(str,v))} for q,bank in banks.items() for b,v in enumerate(bank)])
    for panel in ('state','heldout'):write(f'source_fidelity/{panel}_summary.csv',[r for r in read(ROOT/'source_fidelity/component_statistics.csv') if r['panel']==panel and r['metric']=='source_error'])
    js('source_fidelity/decision.json',{q:json.loads((ROOT/f'source_fidelity/{q}/decision.json').read_text()) for q in SCALES})
def retention(level):
    freeze();r,z,ix,c,g,o,m=arrays();q='q25';samplers={n:Sampler(level,r,z,ix,e,q) for n,e in g.items()};contexts={n:context(s) for n,s in samplers.items()} if level>0 else {};bank=[];audit=[]
    for b in range(B):
        values={}
        for name,s in samplers.items():
            random=rng(f'R{level}|q25|{b}|{name}');attempts=0
            while True:
                v,used=s.draw(random,m);attempts+=1
                if level!=3 or np.array_equal(used,s.targets):break
            values[name]=v
            if level>0:check_incidence(s,contexts[name])
            else:assert np.array_equal(np.bincount(s.last_graph_codes//3704,minlength=3704)[s.nodes],s.size)
            audit.append({'level':level,'draw':b,'orientation':name,'seed':seed(f'R{level}|q25|{b}|{name}'),'conditioning_attempts':attempts,'graph_hash':s.last_graph_hash,'constraints_validated':True})
        bank.append(score(g,o,m,q,values=values))
        if b%256==0:print('R'+str(level),'IID',b,flush=True)
    rel='retention/R'+str(level);obs=score(g,o,m,q);save(rel+'/bank.npz',observed=obs,statistics=np.asarray(bank));write(rel+'/graph_audit.csv',audit);vector_summaries(rel,obs,np.asarray(bank))
    write(f'retention/q25_R{level}.csv',read(ROOT/rel/'component_statistics.csv'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['edges','spatial','matched','degree','retention']);p.add_argument('--variant',default='hash',choices=['hash','pca32','pca64','target_pca32']);p.add_argument('--level',type=int,choices=range(4));a=p.parse_args()
    {'edges':edges,'spatial':spatial,'matched':lambda:matched(a.variant),'degree':degree,'retention':lambda:retention(a.level)}[a.action]()

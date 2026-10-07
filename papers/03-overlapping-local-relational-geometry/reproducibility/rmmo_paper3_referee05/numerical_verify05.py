"""Full read-only replay of all retained draws and normalization authorities."""
import argparse,subprocess,itertools
import scipy.sparse as sp
from common05 import *
from references05 import context,check_incidence
TARGET_CARRIERS={i:set(v) for i,v in enumerate(indices(rows()))}
def close(a,b,label):np.testing.assert_allclose(a,b,rtol=2e-11,atol=2e-12,err_msg=label)
def independent_score(g,o,m,q,codes):
    means={}
    for name,e in g.items():
        nodes=e[q]['nodes'];c=np.asarray(codes[name]);left,right=c//3704,c%3704;assert len(np.unique(c))==len(c)
        assert set(left)==set(nodes) and set(right)<=TARGET_CARRIERS[e['target']]
        counts=np.bincount(np.searchsorted(nodes,left),minlength=len(nodes));assert np.all(counts>0)
        means[name]={}
        incidence=sp.csr_matrix((np.ones(len(c)),(np.searchsorted(nodes,left),right)),shape=(len(nodes),3704))
        for panel,z in m.items():
            s=np.asarray(incidence@z)/counts[:,None];means[name][panel]=unit(s) if panel=='state' else s
    v=[]
    for panel in ('state','heldout'):
        defects=[];errors=[]
        for name,a,b in OVERLAPS:
            nodes=o[name,q];a0=means[a][panel][np.searchsorted(g[a][q]['nodes'],nodes)];b0=means[b][panel][np.searchsorted(g[b][q]['nodes'],nodes)]
            defects.append(np.linalg.norm(a0-b0,axis=1).mean());errors.append(((np.linalg.norm(m[panel][nodes]-a0,axis=1)+np.linalg.norm(m[panel][nodes]-b0,axis=1))/2).mean())
        v+=defects+errors
    return np.array(v)
def fiber_codes(g,q,c=None):
    out={}
    for name,e in g.items():
        if c is None:tar,mem=fiber_matrix(e,q)
        else:
            p='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))));left,right=c[p];tar,mem=fiber_matrix(e,q,right,left) if e['source']<e['target'] else fiber_matrix(e,q,left,right)
        l,t=np.where(mem);out[name]=e[q]['nodes'][l]*3704+tar[t]
    return out
def numerical_family(obs,v,lower):
    # Independent arithmetic, not the report's family routine.
    s=-1 if lower else 1;mu=np.mean(v,axis=0);sd=np.std(v,axis=0,ddof=1);zz=s*(v-mu)/(sd+1e-8);zo=s*(obs-mu)/(sd+1e-8)
    tail=lambda x,y:(1+np.count_nonzero(x>=y))/1025
    comps=[]
    for i in range(len(obs)):
        comps.append({'observed':obs[i],'null_mean':mu[i],'null_sd':sd[i],'absolute_effect':s*(obs[i]-mu[i]),'relative_effect':s*(obs[i]-mu[i])/mu[i] if mu[i] else 0,
            'Z':zo[i],'conditional_tail':tail(s*v[:,i],s*obs[i]),'maxT_tail':tail(zz.max(1),zo[i]),'null95_lo':np.quantile(v[:,i],.025),'null95_hi':np.quantile(v[:,i],.975)})
    return comps,{'T_sum':zo.sum(),'T_min':zo.min(),'sum_tail':tail(zz.sum(1),zo.sum()),'min_tail':tail(zz.min(1),zo.min())}
def check_family(rel,lower=True):
    f=np.load(ROOT/rel/'bank.npz');obs,v=f['observed'],f['statistics'];assert v.shape==(B,len(obs)) and np.isfinite(v).all()
    rows0=read(ROOT/rel/'component_statistics.csv');fr=read(ROOT/rel/'family_statistics.csv');assert len(rows0)==len(obs)
    if len(obs)==3:sets=[(None,None,0)]
    else:sets=[(p,k,j*6+i*3) for j,p in enumerate(('state','heldout')) for i,k in enumerate(('defect','source_error'))]
    for panel,metric,off in sets:
        a,b=numerical_family(obs[off:off+3],v[:,off:off+3],lower);selected=[r for r in rows0 if panel is None or r['panel']==panel and r['metric']==metric]
        for r,w in zip(selected,a):
            for key,val in w.items():close(float(r[key]),val,rel+' '+key)
        fam=next(r for r in fr if panel is None or r['panel']==panel and r['metric']==metric)
        for key,val in b.items():close(float(fam[key]),val,rel+' '+key)
    return obs,v
def incidence_replays(only=None):
    specifications=[('matched_coherence/'+q,'hash','hash',q,2,f'matched|hash|{q}') for q in SCALES]
    specifications +=[(f'representation/{v}','hash' if v=='target_pca32' else v,v,'q25',2,f'matched|{v}|q25') for v in ('pca32','pca64','target_pca32')]
    specifications +=[(f'retention/R{i}','hash','hash','q25',i,f'R{i}|q25') for i in range(4)]
    if only is not None:
        specifications=[s for s in specifications if s[0] in only];assert len(specifications)==len(only)
    for rel,sv,tv,q,level,prefix in specifications:
        r,z,ix,c,g,o,m=arrays(sv,tv);samplers={n:Sampler(level,r,z,ix,e,q) for n,e in g.items()};contexts={n:context(s) for n,s in samplers.items()};obs,v=check_family(rel);close(independent_score(g,o,m,q,fiber_codes(g,q)),obs,rel+' observed')
        audit=read(ROOT/rel/'graph_audit.csv');assert len(audit)==B*6;lookup={(int(r['draw']),r['orientation']):r for r in audit};assert len(lookup)==B*6
        for b in range(B):
            codes={}
            for name,s in samplers.items():
                random=rng(f'{prefix}|{b}|{name}');tries=0
                while True:
                    _,used=s.draw(random,{});tries+=1
                    if level!=3 or np.array_equal(used,s.targets):break
                check_incidence(s,contexts[name]);codes[name]=s.last_graph_codes.copy();ar=lookup[b,name];assert int(ar['seed'])==seed(f'{prefix}|{b}|{name}') and ar['graph_hash']==s.last_graph_hash
                if 'conditioning_attempts' in ar:assert int(ar['conditioning_attempts'])==tries
                assert np.array_equal(np.bincount(codes[name]//3704,minlength=3704)[s.nodes],s.size)
            close(independent_score(g,o,m,q,codes),v[b],rel+' full draw '+str(b))
        print('PASS all 1024 IID incidence graphs and independently scored banks',rel,flush=True)
def degree_replay():
    r,z,ix,c,g,o,m=arrays();draws=np.load(ROOT/'source_fidelity/degree_only_graphs.npz');audit=read(ROOT/'source_fidelity/graph_audit.csv');lookup={(int(r['draw']),r['pair']):r for r in audit};assert len(lookup)==B*3
    banks={q:check_family('source_fidelity/'+q) for q in SCALES}
    for q,(obs,v) in banks.items():close(independent_score(g,o,m,q,fiber_codes(g,q)),obs,'degree observed')
    for b in range(B):
        new={}
        for p,(left,right) in c.items():
            random=rng(f'degree|{p}|{b}');tries=0
            while True:
                t=random.permutation(right);tries+=1
                if len(np.unique(left*3704+t))==len(left):break
            np.testing.assert_array_equal(t,draws[p][b]);np.testing.assert_array_equal(np.sort(t),np.sort(right));assert set(t)==set(right)
            ar=lookup[b,p];assert int(ar['stub_trials'])==tries and int(ar['seed'])==seed(f'degree|{p}|{b}') and ar['graph_hash']==digest_graph(left,t);new[p]=(left,t)
        for q,(obs,v) in banks.items():close(independent_score(g,o,m,q,fiber_codes(g,q,new)),v[b],'degree '+q+' '+str(b))
    print('PASS all 1024 degree-only simple graphs, margins and three-scale independent scores',flush=True)
def spatial_replay():
    r,z,ix,c,g,o,m=arrays();fixed=np.load(ROOT/'spatial_overlap/fixed_cell_strata.npz');dist={i:chord(z[nodes],z[nodes]) for i,nodes in enumerate(ix)}
    for ref,folder in [('S1','S1')]:
        rel='spatial_overlap/'+folder;obs,v=check_family(rel,False);land=np.load(ROOT/rel/'landmark_draws.npz');audit=read(ROOT/rel/'graph_audit.csv');look={(int(r['draw']),r['orientation']):r for r in audit};assert len(look)==B*6
        contexts={}
        for name,e in g.items():
            nodes=ix[e['source']];d=dist[e['source']].copy();np.fill_diagonal(d,np.inf);one=d.min(1);ten=np.partition(d,9,axis=1)[:,9];den=np.searchsorted(np.quantile(ten,[.25,.5,.75]),ten,side='right');nn=np.searchsorted(np.quantile(one,np.arange(.1,1,.1)),one,side='right')
            rad=np.linalg.norm(z[nodes]-z[e['anchors']].mean(0),axis=1);rb=np.searchsorted(np.quantile(rad,np.arange(.1,1,.1)),rad,side='right');keys=[(r[int(n)]['state'],int(den[i]))+((int(rb[i]),int(nn[i])) if ref=='S1' else ()) for i,n in enumerate(nodes)]
            np.testing.assert_array_equal(fixed[ref+'_'+name+'_keys'],[str(k) for k in keys]);pools=defaultdict(list)
            for i,k in enumerate(keys):pools[k].append(i)
            hist=Counter(keys[int(i)] for i in np.searchsorted(nodes,e['anchors']));contexts[name]=(keys,pools,hist)
        for b in range(B):
            dom={}
            for name,e in g.items():
                nodes=ix[e['source']];keys,pools,hist=contexts[name];random=rng(f'spatial|{ref}|{b}|{name}');selected=np.concatenate([random.choice(pools[k],size=n,replace=False) for k,n in hist.items()]);chosen=nodes[selected]
                np.testing.assert_array_equal(chosen,land[name][b]);assert len(set(chosen))==len(chosen) and Counter(keys[int(i)] for i in selected)==hist
                rho=dist[e['source']][:,selected].min(1);radius=np.quantile(rho,.25);dom[name]=nodes[rho<=radius];ar=look[b,name];close(float(ar['radius']),radius,'refit radius');assert int(ar['seed'])==seed(f'spatial|{ref}|{b}|{name}') and ar['landmark_hash']==hashlib.sha256(np.sort(chosen).astype(np.int64).tobytes()).hexdigest()
            np.testing.assert_array_equal([len(np.intersect1d(dom[a],dom[b])) for _,a,b in OVERLAPS],v[b])
        print('PASS all 1024 independently replayed fixed-stratum spatial draws',ref,flush=True)
def edge_replay():
    r,z,ix,c,g,o,m=arrays();_,_,bins=pair_distances(z,ix);states=sorted({x['state'] for x in r});ns=len(states);st=np.array([states.index(x['state']) for x in r]);obs,v=check_family('edge_agreement');ar=read(ROOT/'edge_agreement/graph_audit.csv');lookup={(int(x['draw']),x['pair']):x for x in ar};assert len(lookup)==B*3
    for j,(p,(a0,b0)) in enumerate(ENDS.items()):
        l,t=c[p];keys=(st[l]*ns+st[t])*10+bins[p][np.searchsorted(ix[a0],l),np.searchsorted(ix[b0],t)];poolkeys=(st[ix[a0]][:,None]*ns+st[ix[b0]][None,:])*10+bins[p];ctx=[]
        for key,n in zip(*np.unique(keys,return_counts=True)):
            left,right=np.where(poolkeys==key);ctx.append((int(key),ix[a0][left],ix[b0][right],int(n)))
        close(np.linalg.norm(m['heldout'][l]-m['heldout'][t],axis=1).mean(),obs[j],'edge observed')
        for b in range(B):
            random=rng(f'edge|{p}|{b}');ll=[];tt=[]
            for key,l,t,n in ctx:
                take=random.choice(len(l),size=n,replace=False);ll.extend(l[take]);tt.extend(t[take])
            ll,tt=np.array(ll),np.array(tt);assert len(set(zip(ll,tt)))==len(ll);keys=(st[ll]*ns+st[tt])*10+bins[p][np.searchsorted(ix[a0],ll),np.searchsorted(ix[b0],tt)];assert Counter(keys)==Counter({k:n for k,l,t,n in ctx})
            row=lookup[b,p];assert int(row['seed'])==seed(f'edge|{p}|{b}') and row['graph_hash']==digest_graph(ll,tt);close(np.linalg.norm(m['heldout'][ll]-m['heldout'][tt],axis=1).mean(),v[b,j],'edge null')
    print('PASS all 1024 matched-edge draws per pair and full constraint histograms',flush=True)

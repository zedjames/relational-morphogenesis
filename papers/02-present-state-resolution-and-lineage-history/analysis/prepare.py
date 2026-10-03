"""Outcome-invariant fold caches and source-only conditional strata."""
import itertools
from utils import *

def main():
    bank=Bank(); CACHE.mkdir(parents=True,exist_ok=True)
    contexts=[]
    for held in REPS:
        for fit in [tuple(r for r in REPS if r!=held)]+[(r,) for r in REPS if r!=held]:
            contexts.append({'held':held,'fit':list(fit),'key':'_'.join(fit)+'__'+held})
    state_names=[c for c in CONFIGS if not c.startswith('target_seed')]
    np.save(CACHE/'history.npy',bank.h)
    np.save(CACHE/'rep.npy',bank.rep)
    np.save(CACHE/'s0.npy',bank.s0)
    groups=[]; stratum_rows=[]
    for rep in REPS:
        state=bank.representations([rep])[0]['hash64_seed0']
        for s0 in sorted(set(bank.s0[bank.rep==rep])):
            remaining=np.flatnonzero((bank.rep==rep)&(bank.s0==s0))
            while len(remaining):
                first=remaining[0]; dist=((state[remaining]-state[first])**2).sum(1)
                ix=np.lexsort((remaining,dist))[:32]; g=np.sort(remaining[ix])
                number=len(groups); groups.append(g.tolist())
                for i in g: stratum_rows.append({'index':int(i),'cellHash':bank.rows[i]['cellHash'],'embryo':rep,'S0':s0,'molecularStratum':number,'stratumSize':len(g)})
                remaining=np.setdiff1d(remaining,g,assume_unique=True)
    write('lineage_family_null/conditional_strata.csv',stratum_rows)
    js('lineage_family_null/strata_authority.json',{'groups':groups,'sourceOnly':True,'targetOrHistoryUsed':False,'cells':3704,'movableCells':sum(len(g) for g in groups if len(g)>1),'strata':len(groups),'singletonStrata':sum(len(g)==1 for g in groups),'maximumSize':max(map(len,groups))})
    for ctx in contexts:
        states,targets=bank.representations(ctx['fit']); tr=np.flatnonzero(np.isin(bank.rep,ctx['fit'])); te=np.flatnonzero(bank.rep==ctx['held'])
        d=CACHE/ctx['key']; d.mkdir(exist_ok=True)
        np.save(d/'train.npy',tr); np.save(d/'test.npy',te)
        for name in state_names:
            z=states[name]; configs=['hash64_seed0']+[f'target_seed{i}' for i in range(1,8)] if name=='hash64_seed0' else [name]
            y=np.stack([targets[c if c.startswith('target_seed') else 'canonical'] for c in configs],axis=1)
            e=d/name; e.mkdir(exist_ok=True)
            arrays={'stateTrain':z[tr],'stateTest':z[te],'targetTrain':y[tr],'actual':y[te]}
            mol=distances(z[te],z[tr]); arrays['molecularDistance']=mol
            arrays['molecularOrder']=np.argsort(mol,axis=1,kind='stable').astype(np.int32)
            for key,value in arrays.items(): np.save(e/(key+'.npy'),value)
            (e/'configurations.json').write_text(json.dumps(configs)+'\n')
            if len(ctx['fit'])==2:
                for j,c in enumerate(configs):
                    old=np.load(PREVIOUS/f'inductive_prediction/predictions/{c}_{ctx["held"]}.npz')
                    assert np.max(np.abs(z[tr]-old['stateTrain']))<1e-10
                    assert np.max(np.abs(y[te,j]-old['actual']))<1e-10
                    assert np.max(np.abs(y[tr,j]-old['targetTrain']))<1e-10
        print('Cached invariant fold',ctx['key'],flush=True)
    js('lineage_family_null/cache_spec.json',{'contexts':contexts,'stateGroups':state_names,'float':'float64','fitRepresentationsMatchPrevious':True,'cacheIsDisposable':True})

if __name__=='__main__': main()

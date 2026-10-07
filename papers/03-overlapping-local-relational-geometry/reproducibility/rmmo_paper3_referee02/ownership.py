"""Fixed-domain retained-structure ladder, exact induced-fiber constraints."""
import argparse
from collections import Counter,defaultdict,deque
import hashlib
import numpy as np
from engine import *
from chains import Chain

def arrays():
    rows,z=load_state();ix=indices(rows);core=cores();geo,overlaps=geometry(z,ix,core)
    return rows,z,ix,core,geo,overlaps,{'state':z,'heldout':load_target()}

class Sampler:
    def __init__(self,level,rows,z,ix,e,scale):
        targets,mem=fiber_matrix(e,scale);nodes=e[scale]['nodes'];alltarget=ix[e['target']]
        anchors=np.unique(e['edge_t']);rho=chord(z[alltarget],z[anchors]).min(1)
        dec=np.searchsorted(np.quantile(rho,np.arange(.1,1,.1)),rho,side='right')
        self.key={int(t): (() if level==0 else (rows[int(t)]['state'],) if level==1
                          else (rows[int(t)]['state'],int(d))) for t,d in zip(alltarget,dec)}
        poolnodes=targets if level>=3 else alltarget
        groups=defaultdict(list)
        for t in poolnodes:groups[self.key[int(t)]].append(int(t))
        self.pools={k:np.asarray(v) for k,v in groups.items()}
        batches=defaultdict(list);self.size=mem.sum(1);self.nodes=nodes;self.targets=targets
        self.observed=mem;self.dec=dec;self.alltarget=alltarget;self.carrier_size=len(rows)
        for i,m in enumerate(mem):
            hist=Counter(self.key[int(t)] for t in targets[m])
            for k,n in hist.items():
                assert n<=len(self.pools[k]);batches[k,n].append(i)
        self.batches={k:np.asarray(v) for k,v in batches.items()}
    @staticmethod
    def sample(pool,k,n,random):
        if k==len(pool):return np.broadcast_to(pool,(n,k))
        if k*2>len(pool):return np.asarray([random.choice(pool,size=k,replace=False) for _ in range(n)])
        x=random.integers(len(pool),size=(n,k))
        while True:
            bad=np.any(np.diff(np.sort(x,axis=1),axis=1)==0,axis=1)
            if not bad.any():break
            x[bad]=random.integers(len(pool),size=(bad.sum(),k))
        return pool[x]
    def draw(self,random,matrices):
        result={p:np.zeros((len(self.nodes),m.shape[1])) for p,m in matrices.items()};used=[];codes=[]
        for (key,k),loc in self.batches.items():
            selected=self.sample(self.pools[key],k,len(loc),random);used.append(selected.ravel())
            codes.append((self.nodes[loc,None]*self.carrier_size+selected).ravel())
            for p,m in matrices.items():result[p][loc]+=m[selected].sum(1)
        for p in result:result[p]/=self.size[:,None]
        self.last_graph_codes=np.sort(np.concatenate(codes).astype(np.int64))
        assert len(np.unique(self.last_graph_codes))==len(self.last_graph_codes)
        self.last_graph_hash=hashlib.sha256(self.last_graph_codes.tobytes()).hexdigest()
        return result,np.unique(np.concatenate(used))

def independent_graph_audit():
    """Exact graph identities for the unchanged R0–R3 frozen-seed draws."""
    import shutil
    freeze();rows,z,ix,core,geo,overlaps,matrices=arrays();records=[]
    for level in range(4):
        mobility=[]
        old_path=ROOT/f'ownership_ladder/R{level}_mobility.csv'
        archive=ROOT/f'freeze/amendments/mean_array_mobility/R{level}_mobility.csv'
        if not archive.exists():archive.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(old_path,archive)
        for scale in SCALES:
            for name,e in geo.items():
                sampler=Sampler(level,rows,z,ix,e,scale);hashes=[];differences=[];previous=None
                li,ri=np.where(sampler.observed)
                observed=np.sort(sampler.nodes[li]*len(rows)+sampler.targets[ri])
                for b in range(B):
                    random=rng(f'R{level}|{scale}|{b}|{name}')
                    while True:
                        _,used=sampler.draw(random,{})
                        if level!=3 or np.array_equal(used,sampler.targets):break
                    codes=sampler.last_graph_codes;hashes.append(sampler.last_graph_hash)
                    assert len(codes)==int(sampler.size.sum())
                    np.testing.assert_array_equal(np.bincount(codes//len(rows),minlength=len(rows))[sampler.nodes],sampler.size)
                    difference=0 if previous is None else len(np.setxor1d(codes,previous))
                    if previous is not None:differences.append(difference)
                    previous=codes.copy()
                    records.append({'level':f'R{level}','scale':scale,'orientation':name,'null_index':b,
                        'graph_hash':sampler.last_graph_hash,'edge_count':len(codes),
                        'edge_symmetric_difference_observed':len(np.setxor1d(codes,observed)),
                        'edge_symmetric_difference_previous':difference,'row_degree_mismatch':0,
                        'used_endpoint_count':len(used),'chain_seed':seed(f'R{level}|{scale}|{b}|{name}'),
                        'sampling_authority':'independent conditional finite draw; not a switch chain'})
                counts=Counter(hashes)
                mobility.append({'level':f'R{level}','scale':scale,'orientation':name,'unique_graphs':len(counts),
                    'draws':B,'mobile':len(counts)>1,'authority':'exact frozen-seed fiber graph hashes',
                    'unique_fraction':len(counts)/B,'duplicate_fraction':1-len(counts)/B,
                    'most_frequent_graph_multiplicity':max(counts.values()),
                    'median_adjacent_graph_symmetric_difference':float(np.median(differences))})
            print('Exact independent graph audit',level,scale,flush=True)
        write(f'ownership_ladder/R{level}_mobility.csv',mobility)
    write('ownership_ladder/independent_graph_bank.csv',records)
    js('ownership_ladder/independent_graph_audit.json',{'members_verified':len(records),
        'outcome_banks_changed':False,'seeds_or_constraints_changed':False,
        'mean_array_hashes_superseded_by_exact_graph_hashes':True,
        'reason':'Floating means are not graph identities; direct replay recovers the exact frozen-seed incidence graphs without altering their scientific measurements.'})

def observed_means(geo,scale,matrices,new=None):
    out={}
    for name,e in geo.items():
        if new is None:tar,mem=fiber_matrix(e,scale)
        else:
            p='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))))
            left,right=new[p]
            tar,mem=fiber_matrix(e,scale,right,left) if e['source']<e['target'] else fiber_matrix(e,scale,left,right)
        out[name]={p:means(m,tar,mem,e[scale]['nodes'],e[scale]['nodes']) for p,m in matrices.items()}
    return out

def measure(values,matrices,geo,overlaps,scale,level,index):
    result=[];ori=[]
    for panel,m in matrices.items():
        prepared={name:unit(v[panel]) if panel=='state' else v[panel] for name,v in values.items()}
        for name,e in geo.items():
            err=np.linalg.norm(m[e[scale]['nodes']]-prepared[name],axis=1)
            ori.append({'level':level,'null_index':index,'scale':scale,'panel':panel,'orientation':name,
                'cells':len(err),'mean_error':float(err.mean()),'median_error':float(np.median(err))})
        for o,a,b in OVERLAPS:
            nodes=overlaps[o,scale];la=np.searchsorted(geo[a][scale]['nodes'],nodes);lb=np.searchsorted(geo[b][scale]['nodes'],nodes)
            first,second=prepared[a][la],prepared[b][lb]
            defect=np.linalg.norm(first-second,axis=1)
            ea=np.linalg.norm(m[nodes]-first,axis=1);eb=np.linalg.norm(m[nodes]-second,axis=1);error=(ea+eb)/2
            result.append({'level':level,'null_index':index,'scale':scale,'panel':panel,'overlap':o,'cells':len(nodes),
                'defect':float(defect.mean()),'median_defect':float(np.median(defect)),
                'source_error':float(error.mean()),'median_source_error':float(np.median(error)),
                'first_orientation_mean_error':float(ea.mean()),'second_orientation_mean_error':float(eb.mean())})
    return result,ori

def signature(geo,new,rows,z,ix):
    sig=[]
    for name,e in geo.items():
        targets=ix[e['target']];rho=chord(z[targets],z[np.unique(e['edge_t'])]).min(1)
        dec=np.searchsorted(np.quantile(rho,np.arange(.1,1,.1)),rho,side='right')
        keys=[(rows[int(t)]['state'],int(d)) for t,d in zip(targets,dec)]
        unique={k:i for i,k in enumerate(sorted(set(keys)))}
        cols=np.asarray([unique[k] for k in keys]);indicator=np.eye(len(unique),dtype=np.int32)[cols]
        p='P'+''.join(map(str,sorted((e['source']+1,e['target']+1))))
        left,right=new[p]
        for scale in SCALES:
            tar,mem=fiber_matrix(e,scale,right,left) if e['source']<e['target'] else fiber_matrix(e,scale,left,right)
            sig.append((mem.astype(np.int32)@indicator[np.searchsorted(targets,tar)]).tobytes())
    return b''.join(sig)

def constrained_components(rows,z,ix,core,geo):
    _,_,bins=pair_distances(z,ix);components={};audit=[]
    for p,(a,b) in ENDS.items():
        left,right=core[p];chain=Chain(left,right,rows,bins[p],ix[a],ix[b]);obj=chain.obj
        pairgeo={n:e for n,e in geo.items() if {e['source'],e['target']}=={a,b}}
        expected=signature(pairgeo,{p:core[p]},rows,z,ix)
        start=tuple(right.tolist());seen={start};queue=deque([start]);checked=0
        choices=[(int(i),int(j)) for g in obj.groups for i in g for j in g if i<j]
        while queue:
            state=queue.popleft();current=np.asarray(state);edge=set(zip(left.tolist(),state))
            for i,j in choices:
                si,sj=int(left[i]),int(left[j]);ti,tj=int(current[i]),int(current[j]);checked+=1
                if si==sj or ti==tj or (si,tj) in edge or (sj,ti) in edge:continue
                local_i,local_j=obj.left[i],obj.left[j];new_i=np.searchsorted(ix[b],tj);new_j=np.searchsorted(ix[b],ti)
                if obj.bins[local_i,new_i]!=obj.edge_bins[i] or obj.bins[local_j,new_j]!=obj.edge_bins[j]:continue
                new=current.copy();new[i],new[j]=tj,ti;key=tuple(new.tolist())
                if key in seen:continue
                if signature(pairgeo,{p:(left,new)},rows,z,ix)!=expected:continue
                seen.add(key);queue.append(key)
                assert len(seen)<=10000000,'Frozen state-space ceiling exceeded; requires explicit incomplete enumeration receipt'
        # Slots with repeated source positions may represent the same graph.
        unique={tuple(sorted(zip(left.tolist(),v))):v for v in seen}
        components[p]=np.asarray(list(unique.values()),dtype=right.dtype)
        audit.append({'pair':p,'reachable_slot_states':len(seen),'reachable_graphs':len(unique),
            'candidate_switches_checked':checked,'all_induced_sizes_and_joint_strata_preserved':True,
            'global_admissible_space_enumerated':False,'component_complete':True,
            'authority':'exact reachable constrained component, not all globally feasible assignments'})
        print('R5 exact reachable component',audit[-1],flush=True)
    js('ownership_ladder/R5_component_audit.json',{'pairs':audit,'mobile':all(len(v)>1 for v in components.values()),
        'nesting':'Fiber-incidence R4 and latent landmark R5 are distinct probability spaces.'})
    np.savez_compressed(ROOT/'ownership_ladder/R5_components.npz',**components)
    return components

def run(level):
    freeze();rows,z,ix,core,geo,overlaps,matrices=arrays()
    bank=[];orient=[];mobility=[]
    folder=ROOT/'ownership_ladder';folder.mkdir(parents=True,exist_ok=True)
    if level==6:
        for scale in SCALES:
            rr,oo=measure(observed_means(geo,scale,matrices),matrices,geo,overlaps,scale,'R6',-1);bank+=rr;orient+=oo
    elif level==5:
        comps=constrained_components(rows,z,ix,core,geo)
        selections={p:np.asarray([rng(f'R5|{b}|{p}').integers(len(v)) for b in range(B)]) for p,v in comps.items()}
        for b in range(B):
            new={p:(core[p][0],comps[p][selections[p][b]]) for p in PAIRS}
            for scale in SCALES:
                rr,oo=measure(observed_means(geo,scale,matrices,new),matrices,geo,overlaps,scale,'R5',b);bank+=rr;orient+=oo
        for p,v in comps.items():mobility.append({'level':'R5','scale':'all','orientation':p,'unique_graphs':len(v),
            'draws':B,'mobile':len(v)>1,'authority':'uniform independent draws over exact reachable component'})
    else:
        for scale in SCALES:
            samplers={n:Sampler(level,rows,z,ix,e,scale) for n,e in geo.items()};chains={}
            hashes=defaultdict(list);changes=Counter()
            if level==4:
                for name,e in geo.items():
                    ss=samplers[name];source,target=np.where(ss.observed)
                    # Exact per-source histogram and fiber-incidence column degrees.
                    bins=np.broadcast_to(ss.dec,(len(ix[e['source']]),len(ss.dec))).astype(np.uint8).copy()
                    chains[name]=Chain(ss.nodes[source],ss.targets[target],rows,bins,ix[e['source']],ss.alltarget,group_by_source=False)
            for b in range(B):
                values={}
                for name,ss in samplers.items():
                    random=rng(f'R{level}|{scale}|{b}|{name}')
                    if level<4:
                        while True:
                            out,used=ss.draw(random,matrices)
                            if level!=3 or np.array_equal(used,ss.targets):break
                        values[name]=out
                        # Mean-array hash is diagnostic, not an exact graph hash.
                        hashes[name].append(hashlib.sha256(out['heldout'].tobytes()).hexdigest())
                    else:
                        chain=chains[name];chain.reset();new,diag=chain.step(random,100)
                        nlocal=np.searchsorted(ss.nodes,ix[geo[name]['source']][chain.obj.left])
                        result={}
                        for panel,m in matrices.items():
                            sums=np.zeros((len(ss.nodes),m.shape[1]));np.add.at(sums,nlocal,m[new])
                            result[panel]=sums/ss.size[:,None]
                        values[name]=result;hashes[name].append(diag['graph_hash']);changes[name]+=diag['edge_symmetric_difference_observed']>0
                rr,oo=measure(values,matrices,geo,overlaps,scale,'R'+str(level),b);bank+=rr;orient+=oo
                if b%128==0:print('R'+str(level),scale,b,flush=True)
            for name,h in hashes.items():mobility.append({'level':'R'+str(level),'scale':scale,'orientation':name,
                'unique_graphs':len(set(h)),'draws':B,'mobile':len(set(h))>1,
                'authority':'exact fiber graph hashes' if level==4 else 'mean-array variability, not graph uniqueness'})
    write(f'ownership_ladder/R{level}_nulls_1024.csv' if level<6 else 'ownership_ladder/observed.csv',bank)
    write(f'ownership_ladder/R{level}_orientations.csv',orient)
    if mobility:write(f'ownership_ladder/R{level}_mobility.csv',mobility)
    print('Ownership R'+str(level)+' complete',flush=True)

def summarize():
    obs=read(ROOT/'ownership_ladder/observed.csv');joint=[];component=[];nr=[];summary=[];profiles=[];terminal_rows=[]
    for level in range(6):
        bank=read(ROOT/f'ownership_ladder/R{level}_nulls_1024.csv');lookup={(int(r['null_index']),r['scale'],r['panel'],r['overlap']):r for r in bank}
        mobility=read(ROOT/f'ownership_ladder/R{level}_mobility.csv');mobile=all(r['mobile']=='True' for r in mobility)
        for scale in SCALES:
            exploratory=mobile;reference_tier=0
            if level==4:
                paths=list((ROOT/'mixing_diagnostics').glob('*/ownership_R4_'+scale+'_decision.json'))
                exploratory=False
                if paths:
                    decision=__import__('json').loads(paths[0].read_text());reference_tier=int(decision['terminal']['tier'])
                    exploratory=decision['terminal']['classification']=='PRACTICAL_EXPLORATION_ADEQUATE'
                    stem=str(paths[0])[:-len('_decision.json')]
                    scored=read(stem+f'_tier{reference_tier}_nulls_1024.csv')
                    scores={(int(v['null_index']),int(v['component'])):float(v['statistic']) for v in scored}
                    for b in range(B):
                        for pi,panel in enumerate(('state','heldout')):
                            for oi,(o,_,_) in enumerate(OVERLAPS):
                                row={'level':'R4','scale':scale,'panel':panel,'overlap':o,'null_index':b,
                                    'defect':scores[b,pi*6+oi],'source_error':scores[b,pi*6+3+oi],
                                    'reference_tier':reference_tier,'proposal_sweeps':100*5**reference_tier}
                                lookup[b,scale,panel,o]=row;terminal_rows.append(row)
            for panel in ('state','heldout'):
                for metric in ('defect','source_error'):
                    os=[next(r for r in obs if r['scale']==scale and r['panel']==panel and r['overlap']==o) for o,_,_ in OVERLAPS]
                    null=np.asarray([[float(lookup[b,scale,panel,o][metric]) for o,_,_ in OVERLAPS] for b in range(B)])
                    meta={'level':f'R{level}','scale':scale,'panel':panel,'metric':metric,'mobile':mobile,'exploration_adequate':exploratory,'reference_tier':reference_tier}
                    j,c,n=save_family(f'ownership_R{level}_{scale}_{panel}_{metric}',[float(r[metric]) for r in os],null,[o for o,_,_ in OVERLAPS],metadata=meta)
                    joint.append(j);component+=c;nr+=n
                    summary.append({**meta,'sum_tail':j['sum_tail'],'min_tail':j['min_tail'],
                        'finite_sufficient':exploratory and j['sum_tail']>=.05,
                        'decision':'FINITELY_SUFFICIENT_FOR_THIS_STATISTIC' if exploratory and j['sum_tail']>=.05 else 'IMMOBILE_REFERENCE' if not mobile else 'MIXING_INADEQUATE' if not exploratory else 'NOT_REPRODUCED'})
    write('ownership_ladder/joint_statistics.csv',joint);write('ownership_ladder/individual_statistics.csv',component)
    write('ownership_ladder/synchronized_null_statistics.csv',nr);write('ownership_ladder/retention_response.csv',summary)
    if terminal_rows:write('ownership_ladder/R4_terminal_nulls_1024.csv',terminal_rows)
    for panel in ('state','heldout'):
        for metric in ('defect','source_error'):
            rr=[r for r in summary if r['scale']=='q25' and r['panel']==panel and r['metric']==metric]
            eligible=[r for r in rr if r['exploration_adequate']];sufficient=[r['level'] for r in eligible if r['finite_sufficient']]
            comparable=[(a,b) for a,b in zip(eligible,eligible[1:]) if int(b['level'][1:])<=4 and int(b['level'][1:])-int(a['level'][1:])==1]
            nonmono=any(a['finite_sufficient'] and not b['finite_sufficient'] for a,b in comparable)
            profiles.append({'panel':panel,'metric':metric,'weakest_sufficient_level':sufficient[0] if sufficient else None,
                'sufficient_levels':sufficient,'decision':'NONMONOTONE_RETENTION_RESPONSE' if nonmono else 'FINITE_RETENTION_PROFILE',
                'nonnested_R4_R5':True,'immobile_levels':[r['level'] for r in rr if not r['mobile']],
                'exploration_inadequate_levels':[r['level'] for r in rr if not r['exploration_adequate']],
                'monotonicity_scope':'adjacent comparable R0–R4 levels only; no ordering theorem for latent R5',
                'observed_R6':'deterministic observed object reproduces by construction; not a randomized inferential reference',
                'weakest_is_demonstrated_not_globally_identified':True})
    js('ownership_ladder/minimal_structure_decision.json',{'profiles':profiles,'criterion':'sum-tail finite non-rejection >=0.05 with mobility; not equivalence proof','biological_n':3})
    for panel,filename in [('state','state'),('heldout','heldout')]:
        write(f'ownership_ladder/{filename}_defect_summary.csv',[r for r in component if r['panel']==panel and r['metric']=='defect'])
        write(f'ownership_ladder/{filename}_source_fidelity.csv',[r for r in component if r['panel']==panel and r['metric']=='source_error'])
    plan=__import__('json').loads((ROOT/'freeze/analysis_plan.json').read_text())
    descriptions=dict(plan['ladder']);descriptions['R4']='R3 plus exact global target fiber-incidence degrees; no additional per-target/source-S0 degree margin. See committed conformance amendment.'
    write('ownership_ladder/level_registry.csv',[{'level':k,'retained_structure':v,'members':1 if k=='R6' else B} for k,v in descriptions.items()])
    text_file('ownership_ladder/specification.md','# Ownership ladder\n\nSee the committed plan, its conformance amendments and the final level registry. Fixed source domains and observed overlap cells are used throughout. R0–R3 exact frozen-seed graph identities, row degrees, duplicate counts and adjacent differences are independently replayed; mean-array hashes are not used as graph identities.\n\nR4 fixes induced fiber-incidence degrees, not latent landmark degrees, and does not retain the superseded extra source-S0-conditioned target-degree margin. R4_terminal_nulls_1024.csv and all primary ownership statistics use the actual terminal escalated banks. R4_nulls_1024.csv, R4_orientations.csv and the original R4_mobility.csv retain the100-sweep sensitivity authority; they are not mislabeled as the terminal reference. Exploration-inadequate terminal banks are excluded from sufficiency decisions.\n\nR5 preserves latent degrees AND checks all induced fiber sizes and joint target-state/anchor-bin histograms. Its exact reachable component is not assumed to equal the entire admissible space. The P13 exact marginal is globally complete only for the explicitly declared source-stub-colored space, not a broader histogram-only law. R4/R5 are not one automatically nested probability space. A non-small reproduction tail is finite non-rejection, never an equivalence proof or a unique global owner.')
    print('Ownership profile',profiles,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('level');a=ap.parse_args()
    if a.level=='summarize':summarize()
    elif a.level=='graph-audit':independent_graph_audit()
    else:run(int(a.level))

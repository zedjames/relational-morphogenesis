#!/usr/bin/env python3
"""Fixed-source correspondence and corrected matched coherence banks."""
import argparse
import datetime
from collections import Counter, defaultdict
import numpy as np
from common import (ROOT,B,PAIRS,ENDS,SCALES,OVERLAPS,check_freeze,load_state,
                    load_target,indices,cores,geometry,fiber_matrix,means,score,
                    unit,chord,pair_distances,rng,write,js,sha,npz_parent,family,null_summary)
from rewire import make_rewirer,degree_only

CONVENTIONS = [('state','unit',True),('state','plain',False),('realization','plain',False),('realization','unit',True)]


def record_geometry(prefix,z,ix,core):
    geo,overlaps=geometry(z,ix,core)
    arrays={'state_embedding':z}
    for name,e in geo.items():
        arrays[name+'_anchors']=e['anchors']
        arrays[name+'_edge_source']=e['edge_s'];arrays[name+'_edge_target']=e['edge_t']
        for s in SCALES:
            arrays[name+'_'+s+'_nodes']=e[s]['nodes'];arrays[name+'_'+s+'_active']=e[s]['active']
            arrays[name+'_'+s+'_radius']=np.asarray(e[s]['radius'])
    for (o,s),nodes in overlaps.items():arrays[o+'_'+s+'_nodes']=nodes
    path=ROOT/prefix/'state_geometry.npz';path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,**arrays)
    js(str(prefix)+'/state_geometry_freeze.json',{'geometry_sha256':sha(path),'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'heldout_panel_opened_for_this_geometry':False,'annotations_history_used':False,'counts':{o+'_'+s:len(n) for (o,s),n in overlaps.items()}})
    return geo,overlaps


def compare_geometry(matrices,geo,overlaps,scale,core_override=None):
    fib={}
    for name,e in geo.items():
        if core_override is None: fib[name]=fiber_matrix(e,scale)
        else:
            s,t=e['source'],e['target'];p='P'+''.join(map(str,sorted((s+1,t+1))))
            left,right=core_override[p]
            if s<t: fib[name]=fiber_matrix(e,scale,right,left)
            else: fib[name]=fiber_matrix(e,scale,left,right)
    results=[]
    for o,a,b in OVERLAPS:
        nodes=overlaps[o,scale]
        if not len(nodes):continue
        for panel,conv,normalized in CONVENTIONS:
            m=matrices[panel]
            first=means(m,*fib[a],nodes,geo[a][scale]['nodes'])
            second=means(m,*fib[b],nodes,geo[b][scale]['nodes'])
            results.append({'overlap':o,'scale':scale,'panel':panel,'convention':conv,'cells':len(nodes),**score(first,second,m[nodes],normalized)})
    return results


class MatchedFiber:
    def __init__(self,rows,z,ix,entry,scale,nodes):
        target_ix=ix[entry['target']]
        # All target cells binned by distance to the exact observed target anchors.
        anchors=np.unique(entry['edge_t'])
        rho=chord(z[target_ix],z[anchors]).min(1)
        breaks=np.quantile(rho,np.arange(.1,1.,.1))
        bins=np.searchsorted(breaks,rho,side='right')
        groups=defaultdict(list)
        self.key={}
        for node,decile in zip(target_ix,bins):
            key=(rows[int(node)]['state'],int(decile));groups[key].append(int(node));self.key[int(node)]=key
        self.pools={k:np.asarray(v) for k,v in groups.items()}
        targets,member=fiber_matrix(entry,scale)
        self.nodes=nodes
        self.specs=[]
        balance=[]
        for node,loc in zip(nodes,np.searchsorted(entry[scale]['nodes'],nodes)):
            observed=targets[member[loc]]
            histogram=Counter(self.key[int(t)] for t in observed)
            self.specs.append(list(histogram.items()))
            for key,size in histogram.items():
                assert size<=len(self.pools[key])
                balance.append({'source_node':int(node),'target_S0':key[0],'target_anchor_decile':key[1],'observed_count':size,'pool_size':len(self.pools[key]),'null_count_per_draw':size,'variable_membership_possible':size<len(self.pools[key])})
        self.balance=balance

    def draw(self,random,matrices):
        output={p:np.zeros((len(self.nodes),m.shape[1])) for p,m in matrices.items()}
        for i,spec in enumerate(self.specs):
            chosen=np.concatenate([random.choice(self.pools[k],size=n,replace=False) for k,n in spec])
            assert len(np.unique(chosen))==len(chosen)
            assert Counter(self.key[int(t)] for t in chosen)==dict(spec)
            for p,m in matrices.items():output[p][i]=m[chosen].mean(0)
        return output


def summarize(prefix,observed,banks):
    summary=[];families=[]
    for name,bank in banks.items():
        lookup={(r['null_index'],r['scale'],r['overlap'],r['panel'],r['convention']):r for r in bank}
        for panel,conv,_ in CONVENTIONS:
            for s in sorted({r['scale'] for r in bank}):
                os=[r for r in observed if r['scale']==s and r['panel']==panel and r['convention']==conv]
                if len(os)!=3:continue
                for metric in ('defect','source_error'):
                    null=np.asarray([[lookup[b,s,o['overlap'],panel,conv][metric] for o in os] for b in range(B)])
                    result=family([o[metric] for o in os],null,lower=True)
                    families.append({'bank':name,'scale':s,'panel':panel,'convention':conv,'metric':metric,**{k:v for k,v in result.items() if k!='maxT_adjusted_tails'},'individual_maxT_tails':','.join(map(str,result['maxT_adjusted_tails']))})
                    for j,o in enumerate(os):
                        summary.append({'bank':name,'overlap':o['overlap'],'scale':s,'panel':panel,'convention':conv,'metric':metric,'cells':o['cells'],**null_summary(o[metric],null[:,j],lower=True),'maxT_tail':result['maxT_adjusted_tails'][j]})
    write(str(prefix)+'/summary.csv',summary)
    write(str(prefix)+'/familywise_summary.csv',families)
    return summary,families


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--family',default='hash_family_reference');ap.add_argument('--matched-only',action='store_true');ap.add_argument('--target-pca',action='store_true');args=ap.parse_args()
    if args.target_pca:args.family='target_pca_control'
    check_freeze();rows,z=load_state();ix=indices(rows);core=cores()
    prefix='representation_families/'+args.family
    if args.family not in ('hash_family_reference','target_pca_control'):
        with np.load(ROOT/prefix/'representation.npz') as f:z=f['embedding']
        core={}
        for p,(a,b) in ENDS.items():
            block=chord(z[ix[a]],z[ix[b]])
            epsilon=np.quantile(np.r_[block.min(1),block.min(0)],.5)
            left,right=np.where(block<=epsilon);core[p]=(ix[a][left],ix[b][right])
        js(prefix+'/pair_calibration.json',{p:{'edges':len(c[0]),'source_landmarks':len(np.unique(c[0])),'target_landmarks':len(np.unique(c[1])),'stable_core_7_of_9_applied':False} for p,c in core.items()})
    geo,overlaps=record_geometry(prefix,z,ix,core)
    if args.family=='hash_family_reference':
        old=npz_parent(16,'overlaps/overlap_domains.npz')
        for (o,s),nodes in overlaps.items():assert np.array_equal(nodes,old[o+'_'+s]),(o,s)
    # ZR opened only after state geometry has a persisted hash receipt.
    if args.target_pca:
        with np.load(ROOT/prefix/'representation.npz') as f:target=f['embedding']
    else:target=load_target()
    matrices={'state':z,'realization':target}
    observed=[r for scale in SCALES for r in compare_geometry(matrices,geo,overlaps,scale)]
    write(prefix+'/observed_coherence.csv',observed)
    js(prefix+'/target_reveal_receipt.json',{'state_geometry_sha256_verified':sha(ROOT/prefix/'state_geometry.npz'),'target_reveal_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'panel_term':'held-out contemporaneous molecular panel','biological_n':3})
    matchers={};balance=[]
    for o,a,b in OVERLAPS:
        for s in SCALES:
            nodes=overlaps[o,s]
            for name in (a,b):
                mm=MatchedFiber(rows,z,ix,geo[name],s,nodes);matchers[o,s,name]=mm
                balance.extend({'overlap':o,'scale':s,'orientation':name,**r} for r in mm.balance)
    write(prefix+'/matched_balance.csv',balance)
    matched=[]
    for b in range(B):
        for o,a,c in OVERLAPS:
            for s in SCALES:
                nodes=overlaps[o,s]
                if not len(nodes):continue
                left=matchers[o,s,a].draw(rng(f'matched|{args.family}|{b}|{o}|{s}|{a}'),matrices)
                right=matchers[o,s,c].draw(rng(f'matched|{args.family}|{b}|{o}|{s}|{c}'),matrices)
                for panel,conv,norm in CONVENTIONS:
                    matched.append({'null_index':b,'overlap':o,'scale':s,'panel':panel,'convention':conv,'cells':len(nodes),**score(left[panel],right[panel],matrices[panel][nodes],norm),'stratum_mismatch':0,'fiber_size_mismatch':0})
        if b%128==0:print(args.family,'matched',b,flush=True)
    write(prefix+'/matched_nulls_1024.csv',matched)
    banks={'corrected_matched':matched}
    if not args.matched_only:
        _,_,bins=pair_distances(z,ix)
        objects={p:make_rewirer(*core[p],rows,bins[p],ix[a],ix[c]) for p,(a,c) in ENDS.items()}
        for kind in ('strict_S0_distance_degree','degree_only_uniform_stub_rejection'):
            bank=[];mobility=[];edges={p:[] for p in PAIRS}
            for b in range(B):
                new={}
                for p in PAIRS:
                    left,right=core[p]
                    random=rng(f'correspondence|{args.family}|{kind}|{b}|{p}')
                    if kind.startswith('strict'):target,diag=objects[p].draw(random,sweeps=10)
                    else:target,diag=degree_only(left,right,random)
                    new[p]=(left,target);edges[p].append(target)
                    mobility.append({'null_index':b,'pair':p,'bank':kind,**diag})
                bank.extend({'null_index':b,**r} for r in compare_geometry(matrices,geo,overlaps,'q25',new))
                if b%128==0:print(args.family,kind,b,flush=True)
            write(prefix+'/'+kind+'_nulls_1024.csv',bank)
            write(prefix+'/'+kind+'_mobility.csv',mobility)
            np.savez_compressed(ROOT/prefix/(kind+'_rewired_cores.npz'),**{p+'_target':np.asarray(v) for p,v in edges.items()})
            banks[kind]=bank
    summaries,families=summarize(prefix,observed,banks)
    js(prefix+'/decision.json',{'banks':list(banks),'replicates_per_bank':B,'biological_n':3,'all_three_family_results':families,'interpretation':'conditional finite-carrier sensitivity; strict-bank mobility is required for correspondence-breaking authority'})
    print(args.family,'complete',flush=True)


if __name__=='__main__':main()

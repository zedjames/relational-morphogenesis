"""Stage B bindings: reused immutable numerics, own inputs/seeds/output directory."""
import os,sys
from collections import Counter,defaultdict
from io05 import *
P1=ROOT.parent/'rmmo_paper3_referee01';P2=ROOT.parent/'rmmo_paper3_referee02';P3=ROOT.parent/'rmmo_paper3_referee03';P4=ROOT.parent/'rmmo_paper3_referee04'
sys.path.insert(0,str(P2));sys.path.insert(0,str(P1))
import common as historical
from ownership import Sampler,measure,observed_means
indices,geometry,fiber_matrix,pair_distances=historical.indices,historical.geometry,historical.fiber_matrix,historical.pair_distances
def revealed():
    r=json.loads((ROOT/'separation_audit/reveal_receipt.json').read_text());assert r['loaded_only_after_geometry_commit'];assert sha(ROOT/'geometry/stageA_arrays.npz')==r['geometry_sha256'];return r
def rows(full=False):
    revealed();r=read(ROOT/'authority/cell_carrier.csv.gz');return r if full else [{k:x[k] for k in ('cellHash','replicate','state')} for x in r]
def state(variant='hash'):
    if variant in ('hash','hash_primary'):return np.load(ROOT/'strict_lanes/lane_embeddings/S0.npz')['embedding']
    revealed();return unit(np.load(ROOT/f'evaluation_inputs/{variant}/representation.npz')['embedding'])
def target(variant='hash'):
    revealed();name='target_pca32' if variant=='target_pca32' else 'heldout_hash';return np.load(ROOT/f'evaluation_inputs/{name}/representation.npz')['embedding'].astype(float)
def arrays(variant='hash',panel='hash'):
    r=rows();z=state(variant);ix=indices(r);path='geometry/pair_cores.npz' if variant=='hash' else f'geometry/{variant}/pair_cores.npz';f=np.load(ROOT/path);c={p:(f[p+'_left'],f[p+'_right']) for p in PAIRS};g,o=geometry(z,ix,c)
    return r,z,ix,c,g,o,{'state':z,'heldout':target(panel)}
def score(g,o,m,q,edges=None,values=None):
    rr,_=measure(values if values is not None else observed_means(g,q,m,edges),m,g,o,q,'direct',0);return np.array([r[metric] for panel in ('state','heldout') for metric in ('defect','source_error') for r in rr if r['panel']==panel])
def digest_graph(left,right):return hashlib.sha256(np.sort(np.asarray(left,dtype=np.int64)*3704+right).tobytes()).hexdigest()
def family(obs,v,lower=True):
    v=np.asarray(v,float);obs=np.asarray(obs,float);assert v.shape==(B,len(obs)) and np.isfinite(v).all();s=-1 if lower else 1;mu=v.mean(0);sd=v.std(0,ddof=1);z=s*(v-mu)/(sd+1e-8);zo=s*(obs-mu)/(sd+1e-8);tail=lambda x,y:float((1+np.count_nonzero(x>=y))/(B+1))
    rr=[{'observed':float(x),'null_mean':float(mu[i]),'null_sd':float(sd[i]),'null95_lo':float(np.quantile(v[:,i],.025)),'null95_hi':float(np.quantile(v[:,i],.975)),'absolute_effect':float(s*(x-mu[i])),'relative_effect':float(s*(x-mu[i])/mu[i]) if mu[i] else 0.,'Z':float(zo[i]),'conditional_tail':tail(s*v[:,i],s*x),'maxT_tail':tail(z.max(1),zo[i])} for i,x in enumerate(obs)]
    return {'T_sum':float(zo.sum()),'T_min':float(zo.min()),'sum_tail':tail(z.sum(1),zo.sum()),'min_tail':tail(z.min(1),zo.min()),'componentwise_family_supported':all(r['Z']>0 and r['maxT_tail']<=.05 for r in rr),'components':rr}
def summarize(rel,obs,bank,labels,lower=True,meta=None):
    meta=meta or {};d=family(obs,bank,lower);write(rel+'/component_statistics.csv',[{**meta,'component':label,**r,'T_sum':d['T_sum'],'T_min':d['T_min'],'sum_tail':d['sum_tail'],'min_tail':d['min_tail']} for label,r in zip(labels,d['components'])]);write(rel+'/family_statistics.csv',[{**meta,**{k:v for k,v in d.items() if k!='components'}}]);js(rel+'/decision.json',d);return d
def vector_summaries(rel,obs,bank):
    comps=[];fams=[];dec={}
    for j,panel in enumerate(('state','heldout')):
        for k,metric in enumerate(('defect','source_error')):
            off=j*6+k*3;d=family(obs[off:off+3],bank[:,off:off+3]);dec[panel+'_'+metric]=d;comps.extend({'panel':panel,'metric':metric,'component':o,**r,'T_sum':d['T_sum'],'T_min':d['T_min'],'sum_tail':d['sum_tail'],'min_tail':d['min_tail']} for (o,_,_),r in zip(OVERLAPS,d['components']));fams.append({'panel':panel,'metric':metric,**{k:v for k,v in d.items() if k!='components'}})
    write(rel+'/component_statistics.csv',comps);write(rel+'/family_statistics.csv',fams);js(rel+'/decision.json',dec);return dec

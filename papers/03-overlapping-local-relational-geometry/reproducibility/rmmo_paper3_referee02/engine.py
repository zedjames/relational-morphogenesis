"""Owned outputs; immutable parent numerical functions only, never parent writers."""
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent / 'rmmo_paper3_referee01'
sys.path.insert(0, str(PARENT))
import common as old

B = 1024
PAIRS, ENDS, OVERLAPS, SCALES = old.PAIRS, old.ENDS, old.OVERLAPS, old.SCALES
read, sha, unit, chord = old.read, old.sha, old.unit, old.chord
load_state, load_target, indices, cores = old.load_state, old.load_target, old.indices, old.cores
geometry, fiber_matrix, means, pair_distances = old.geometry, old.fiber_matrix, old.means, old.pair_distances
npz_parent = old.npz_parent

def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def seed(label): return int.from_bytes(hashlib.sha256(('RMMO_P3_REFEREE02_V1|'+label).encode()).digest()[:8], 'big')
def rng(label): return np.random.default_rng(seed(label))
def write(rel, rows, fields=None):
    rows=list(rows); path=ROOT/rel; path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or (list(rows[0]) if rows else ['status'])
    with (gzip.open if path.suffix=='.gz' else open)(path,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def js(rel, obj):
    p=ROOT/rel;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
def text_file(rel, value):
    p=ROOT/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(value.rstrip()+'\n')
def freeze():
    r=json.loads((ROOT/'freeze/freeze_receipt.json').read_text())
    for p,h in r['file_sha256'].items():assert sha(ROOT/p)==h,('plan changed',p)
    return r
def display(p): return '<=0.001' if abs(p-1/1025)<1e-12 else f'{p:.3f}'
def family(obs, null, lower=True):
    v=np.asarray(null,float);o=np.asarray(obs,float)
    assert v.shape==(B,len(o)) and np.isfinite(v).all() and np.isfinite(o).all()
    sign=-1 if lower else 1;mu=v.mean(0);sd=v.std(0,ddof=1)
    z=sign*(v-mu)/(sd+1e-8);zo=sign*(o-mu)/(sd+1e-8)
    tail=lambda vv,oo:float((1+np.sum(vv>=oo))/(B+1))
    component=[]
    for j,x in enumerate(o):
        p=tail(sign*v[:,j],sign*x);adj=tail(z.max(1),zo[j])
        component.append({'observed':float(x),'null_mean':float(mu[j]),'null_sd':float(sd[j]),
            'null95':float(np.quantile(v[:,j],.95)),'Z':float(zo[j]),'conditional_tail':p,
            'maxT_tail':adj,'tail_display':display(p),'absolute_effect':float(mu[j]-x),
            'relative_effect':float((mu[j]-x)/mu[j]) if mu[j] else 0.})
    ps=tail(z.sum(1),zo.sum());pm=tail(z.min(1),zo.min())
    return {'T_sum':float(zo.sum()),'T_min':float(zo.min()),'sum_tail':ps,'min_tail':pm,
        'sum_tail_display':display(ps),'min_tail_display':display(pm),
        'joint_aggregate_supported':ps<=.05,'weakest_component_supported':pm<=.05,
        'componentwise_family_supported':all(x['maxT_tail']<=.05 and x['Z']>0 for x in component),
        'components':component},z

def save_family(name,obs,null,labels,lower=True,metadata=None):
    result,z=family(obs,null,lower);meta=metadata or {}
    joint={**meta,'family':name,**{k:v for k,v in result.items() if k!='components'}}
    components=[{**meta,'family':name,'component':label,**r,'T_sum':result['T_sum'],
        'T_min':result['T_min'],'sum_tail':result['sum_tail'],'min_tail':result['min_tail']}
        for label,r in zip(labels,result['components'])]
    nullrows=[{**meta,'family':name,'null_index':b,'T_sum':float(z[b].sum()),'T_min':float(z[b].min()),
              'raw_components':'|'.join(map(str,null[b])), 'standardized_components':'|'.join(map(str,z[b]))}
              for b in range(B)]
    return joint,components,nullrows

"""Authorized calibration-only banks, frozen before output; immutable inferential inputs."""
import argparse,sys
from collections import Counter
from io08 import *
P5=ROOT.parent/'rmmo_paper3_referee05';P7=ROOT.parent/'rmmo_paper3_acceptance07';sys.path.insert(0,str(P5));import common05 as p
NAMESPACE='RMMO_P3_FINAL08_CALIBRATION_V1|'
def seed(label):return int.from_bytes(hashlib.sha256((NAMESPACE+label).encode()).digest()[:8],'big')
def rng(label):return np.random.default_rng(seed(label))
def freeze():
    plan=json.loads((ROOT/'freeze/missing_calibration_plan.json').read_text());rec=json.loads((ROOT/'freeze/missing_calibration_commit.json').read_text());assert sha(ROOT/'freeze/missing_calibration_plan.json')==rec['plan_sha256'];assert sha(Path(__file__))==rec['implementation_sha256'];return plan
class Sampler(p.Sampler):
    @staticmethod
    def sample(pool,k,n,random):
        if k>8 and k!=len(pool):return np.asarray([random.choice(pool,size=k,replace=False) for _ in range(n)])
        return p.Sampler.sample(pool,k,n,random)
def context(s):
    keys={k:i for i,k in enumerate(sorted(set(s.key.values())))};codes=s.nodes[np.where(s.observed)[0]]*3704+s.targets[np.where(s.observed)[1]];lookup=np.array([keys[s.key[t]] if t in s.key else -1 for t in range(3704)]);return len(keys),lookup,np.sort(codes//3704*len(keys)+lookup[codes%3704])
def check(s,ctx):
    n,key,expected=ctx;c=s.last_graph_codes;np.testing.assert_array_equal(np.sort(c//3704*n+key[c%3704]),expected);assert len(np.unique(c))==len(c);np.testing.assert_array_equal(np.bincount(c//3704,minlength=3704)[s.nodes],s.size)
def save(rel,**a):
    path=ROOT/'independent_calibration_bank'/rel;path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**a)
def inferential(bank):return (P7 if bank.startswith('acceptance07/') else P5)/bank/'bank.npz'
def vector(g,o,m,q,metric,values=None,edges=None):
    off=0 if metric=='defect' else 3;return p.score(g,o,m,q,values=values,edges=edges)[[off,off+1,off+2,off+6,off+7,off+8]]
def obs(bank,g,o,m,q,metric):
    f=np.load(inferential(bank));off=0 if metric=='defect' else 3;expected=f['observed'] if len(f['observed'])==6 else f['observed'][[off,off+1,off+2,off+6,off+7,off+8]];np.testing.assert_allclose(vector(g,o,m,q,metric),expected,rtol=2e-11,atol=2e-12);return expected
def matched(bank):
    plan=freeze();spec=next(x for x in plan['banks'] if x['bank']==bank);assert spec['metric']=='defect';assert not (ROOT/'independent_calibration_bank'/bank/'bank.npz').exists(),'No overwrite, rescue or bank selection';variant,q=spec['variant'],spec['scale'];r,z,ix,c,g,o,m=p.arrays(variant);samplers={n:Sampler(2,r,z,ix,e,q) for n,e in g.items()};ctx={n:context(s) for n,s in samplers.items()};v=[];audit=[];observed=obs(bank,g,o,m,q,'defect')
    for b in range(B):
        values={}
        for name,s in samplers.items():
            label=f'matched|{bank}|{b}|{name}';values[name],_=s.draw(rng(label),m);check(s,ctx[name]);audit.append(dict(draw=b,orientation=name,seed=seed(label),graph_hash=s.last_graph_hash,matching_validated=True,fallback=0,replacement=False,NOT_INFERENTIAL=True))
        v.append(vector(g,o,m,q,'defect',values=values))
        if b%256==0:print('Calibration-only matched draw',bank,b,'of',B,flush=True)
    save(bank+'/bank.npz',observed=observed,statistics=np.asarray(v));write('independent_calibration_bank/'+bank+'/graph_audit.csv',audit);print('Calibration-only bank saved',bank,'no inferential rank inspected',flush=True)
def degree():
    plan=freeze();specs=[x for x in plan['banks'] if x['metric']=='source_error'];assert len(specs)==2
    for x in specs:assert not (ROOT/'independent_calibration_bank'/x['bank']/'bank.npz').exists()
    r,z,ix,c,g,o,m=p.arrays();v={x['scale']:[] for x in specs};edges={name:[] for name in PAIRS};audit=[];observed={x['scale']:obs(x['bank'],g,o,m,x['scale'],'source_error') for x in specs}
    for b in range(B):
        new={}
        for name,(left,right) in c.items():
            label=f'degree_missing|{name}|{b}';random=rng(label);tries=0
            while True:
                t=random.permutation(right);tries+=1
                if len(np.unique(left*3704+t))==len(left):break
            np.testing.assert_array_equal(np.sort(t),np.sort(right));new[name]=(left,t);edges[name].append(t);audit.append(dict(draw=b,pair=name,seed=seed(label),stub_trials=tries,graph_hash=p.digest_graph(left,t),simplicity_degrees_endpoints_validated=True,NOT_INFERENTIAL=True))
        for q in v:v[q].append(vector(g,o,m,q,'source_error',edges=new))
        if b%256==0:print('Calibration-only shared degree draw',b,'of',B,flush=True)
    save('source_fidelity/degree_only_graphs.npz',**{name:np.asarray(vv) for name,vv in edges.items()});write('independent_calibration_bank/source_fidelity/graph_audit.csv',audit)
    for x in specs:save(x['bank']+'/bank.npz',observed=observed[x['scale']],statistics=np.asarray(v[x['scale']]))
    print('Two degree scale-calibration banks saved; no inferential rank inspected',flush=True)
if __name__=='__main__':
    args=argparse.ArgumentParser();args.add_argument('action',choices=['matched','degree']);args.add_argument('--bank');a=args.parse_args();matched(a.bank) if a.action=='matched' else degree()

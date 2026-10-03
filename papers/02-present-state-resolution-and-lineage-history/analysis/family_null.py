"""Exact nested replay of the complete prospective lineage family.

The C++ nearest-neighbor kernel prunes only certified nonnegative lower bounds.
Ridge uses a seven-feature Schur complement of the exact normal equations.
Immutable outcome-invariant fold caches may be regenerated with prepare.py.
"""
import argparse
import ctypes
import time
from scipy.linalg import cho_factor, cho_solve
from utils import *

STAGES=['original','E1']
MODELS=['jointMetric','ridge']

def nse_many(actual,pred):
    a=actual/np.where(np.linalg.norm(actual,axis=2,keepdims=True)==0,1,np.linalg.norm(actual,axis=2,keepdims=True))
    p=pred/np.where(np.linalg.norm(pred,axis=2,keepdims=True)==0,1,np.linalg.norm(pred,axis=2,keepdims=True))
    return np.mean(((a-p)**2).sum(2)/np.maximum((a*a).sum(2),1e-12),axis=0)

class Kernel:
    def __init__(self):
        self.lib=ctypes.CDLL(str(CACHE/'nearest.so'))
        i=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        f=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        self.lib.metric_scores.argtypes=[ctypes.c_int]*4+[i,f,f,i,i,f,f,ctypes.c_double,f]
        self.lib.metric_scores.restype=None
    def run(self,ctx,hd,lam):
        n,nt=ctx.mol.shape; q,d=ctx.actual.shape[1:]
        out=np.empty((q,3)); pi=np.arange(n,dtype=np.int32); pj=np.arange(nt,dtype=np.int32)
        self.lib.metric_scores(n,nt,q,d,ctx.order,ctx.mol,np.ascontiguousarray(hd),pi,pj,ctx.y,ctx.actual,lam,out)
        assert np.all(np.isfinite(out))
        return out

class Context:
    def __init__(self,spec,name):
        self.spec=spec; self.name=name; p=CACHE/spec['key']; e=p/name
        self.tr=np.load(p/'train.npy'); self.te=np.load(p/'test.npy')
        self.configs=json.loads((e/'configurations.json').read_text())
        get=lambda s:np.load(e/(s+'.npy'),mmap_mode='r')
        self.x,self.xt=standard(get('stateTrain'),get('stateTest'))
        self.y=get('targetTrain'); self.actual=get('actual'); self.mol=get('molecularDistance'); self.order=get('molecularOrder')
        self.mean=self.y.mean(0); yy=self.y.reshape(len(self.tr),-1)-self.mean.ravel()
        self.q=self.y.shape[1]; self.dim=self.y.shape[2]
        self.ch=[]; self.coef=[]; self.pred=[]; self.ridgeNSE=[]
        xx=self.x.T@self.x; xy=self.x.T@yy
        for alpha in ALPHA:
            ch=cho_factor(xx+alpha*np.eye(xx.shape[0])); coef=cho_solve(ch,xy)
            pred=(self.mean.ravel()+self.xt@coef).reshape(self.actual.shape)
            self.ch.append(ch); self.coef.append(coef); self.pred.append(pred)
            self.ridgeNSE.append(nse_many(self.actual,pred))
        self.ridgeNSE=np.array(self.ridgeNSE).T
    def ridge_history(self,h):
        a,b=standard(history(h[self.tr]),history(h[self.te]))
        xh=self.x.T@a; hh=a.T@a
        hy=a.T@(self.y.reshape(len(self.tr),-1)-self.mean.ravel())
        out=[]
        for j,alpha in enumerate(ALPHA):
            v=cho_solve(self.ch[j],xh)
            bh=np.linalg.solve(hh+alpha*np.eye(7)-xh.T@v,hy-xh.T@self.coef[j])
            pred=self.pred[j].reshape(len(self.te),-1)+(b-self.xt@v)@bh
            out.append(nse_many(self.actual,pred.reshape(self.actual.shape)))
        return np.array(out).T

def choose(inner,allowed,kind):
    if kind=='metric':
        return min([(k,l) for k in range(3) for l in allowed],key=lambda p:(float(inner[p]),K[p[0]],LAM[p[1]]))
    return (min(allowed,key=lambda a:(float(inner[a]),ALPHA[a])),)

class Replay:
    def __init__(self):
        self.spec=json.loads((ROOT/'lineage_family_null/cache_spec.json').read_text())
        self.contexts=[Context(c,name) for c in self.spec['contexts'] for name in self.spec['stateGroups']]
        self.kernel=Kernel(); self.h=np.load(CACHE/'history.npy')
        self.groups=json.loads((ROOT/'lineage_family_null/strata_authority.json').read_text())['groups']
        self.config_index={c:i for i,c in enumerate(CONFIGS)}
        self.members=[(s,c,m) for s in STAGES for c in CONFIGS for m in MODELS]
        self.base_inner_m=np.empty((6,19,3)); self.base_inner_r=np.empty((6,19,7))
        self.inner_spec=[c for c in self.spec['contexts'] if len(c['fit'])==1]
        self.inner_index={c['key']:i for i,c in enumerate(self.inner_spec)}
        self.outer={}; self.baseline={}
        for ctx in self.contexts:
            if len(ctx.spec['fit'])==1:
                b=self.kernel.run(ctx,np.zeros_like(ctx.mol),0); row=self.inner_index[ctx.spec['key']]
                for j,c in enumerate(ctx.configs):
                    ix=self.config_index[c]; self.base_inner_m[row,ix]=b[j]; self.base_inner_r[row,ix]=ctx.ridgeNSE[j]
            else: self.outer[(ctx.spec['held'],ctx.name)]=ctx
        self.mol_params={}
        for e,held in enumerate(REPS):
            rows=self.inner_rows(held)
            for ci,c in enumerate(CONFIGS):
                name='hash64_seed0' if c.startswith('target_seed') else c; ctx=self.outer[(held,name)]; j=ctx.configs.index(c)
                k=min(range(3),key=lambda k:(self.base_inner_m[rows,ci,k].mean(),K[k]))
                base_m=self.kernel.run(ctx,np.zeros_like(ctx.mol),0)[j,k]
                for si,stage in enumerate(STAGES):
                    aa=list(range(4)) if stage=='original' else list(range(7))
                    alpha=min(aa,key=lambda a:(self.base_inner_r[rows,ci,a].mean(),ALPHA[a]))
                    self.baseline[(si,ci,e,0)]=base_m
                    self.baseline[(si,ci,e,1)]=ctx.ridgeNSE[j,alpha]
                    self.mol_params[(si,ci,e)]=(k,alpha)
    def inner_rows(self,held):
        return [i for i,c in enumerate(self.inner_spec) if held!=c['held'] and held not in c['fit']]
    def permutation(self,b):
        ix=np.arange(len(self.h)); rng=np.random.default_rng(stable(f'RMMO_P2_REFEREE02_LINEAGE_NULL_{b}'))
        for g in self.groups: ix[g]=rng.permutation(g)
        return ix
    def run(self,b):
        ix=np.arange(len(self.h)) if b<0 else self.permutation(b); h=self.h[ix]
        im=np.empty((6,19,3,7)); ir=np.empty((6,19,7)); hd_cache={}
        for ctx in self.contexts:
            if len(ctx.spec['fit'])!=1: continue
            key=ctx.spec['key']; row=self.inner_index[key]
            if key not in hd_cache:
                a,v=standard(history(h[ctx.tr]),history(h[ctx.te])); hd_cache[key]=distances(v,a)/7
            hd=hd_cache[key]; mm=np.array([self.kernel.run(ctx,hd,l) for l in LAM]).transpose(1,2,0)
            rr=ctx.ridge_history(h)
            for j,c in enumerate(ctx.configs):
                ci=self.config_index[c]; im[row,ci]=mm[j]; ir[row,ci]=rr[j]
        params=np.empty((2,19,3,3),dtype=np.int32)
        gains=np.empty((2,19,2,3)); outer_scores=np.empty_like(gains)
        for e,held in enumerate(REPS):
            rows=self.inner_rows(held)
            for ci,c in enumerate(CONFIGS):
                for si,stage in enumerate(STAGES):
                    ll=[0,3,4,5] if stage=='original' else list(range(7)); aa=list(range(4)) if stage=='original' else list(range(7))
                    km=choose(im[rows,ci].mean(0),ll,'metric'); ar=choose(ir[rows,ci].mean(0),aa,'ridge')[0]
                    params[si,ci,e]=[*km,ar]
            for name in self.spec['stateGroups']:
                ctx=self.outer[(held,name)]; key=ctx.spec['key']
                if key not in hd_cache:
                    a,v=standard(history(h[ctx.tr]),history(h[ctx.te])); hd_cache[key]=distances(v,a)/7
                ll=sorted({int(params[si,self.config_index[c],e,1]) for c in ctx.configs for si in range(2)})
                mm={li:self.kernel.run(ctx,hd_cache[key],LAM[li]) for li in ll}
                rr=ctx.ridge_history(h)
                for j,c in enumerate(ctx.configs):
                    ci=self.config_index[c]
                    for si in range(2):
                        ki,li,ai=params[si,ci,e]
                        for mi,value in enumerate([mm[int(li)][j,int(ki)],rr[j,int(ai)]]):
                            outer_scores[si,ci,mi,e]=value
                            gains[si,ci,mi,e]=self.baseline[(si,ci,e,mi)]-value
        assert np.all(np.isfinite(gains)) and np.all(np.isfinite(im)) and np.all(np.isfinite(ir))
        out=ROOT/'lineage_family_null/replay'; out.mkdir(exist_ok=True)
        final=out/('observed.npz' if b<0 else f'null_{b:04d}.npz')
        pending=out/('pending_observed.npz' if b<0 else f'pending_{b:04d}.npz')
        np.savez_compressed(pending,assignment=ix,innerMetricNSE=im,innerRidgeNSE=ir,selected=params,gains=gains,outerHistoryNSE=outer_scores)
        os.replace(pending,final)
        return gains

def main():
    p=argparse.ArgumentParser(); p.add_argument('--start',type=int,default=-1); p.add_argument('--stop',type=int,default=999); p.add_argument('--step',type=int,default=1); p.add_argument('--resume',action='store_true'); args=p.parse_args()
    t=time.monotonic(); replay=Replay(); print('Initialized exact replay in',round(time.monotonic()-t,2),'seconds',flush=True)
    for b in range(args.start,args.stop,args.step):
        path=ROOT/'lineage_family_null/replay'/('observed.npz' if b<0 else f'null_{b:04d}.npz')
        if args.resume and path.exists(): continue
        t=time.monotonic(); g=replay.run(b)
        print('Replay',b,'seconds',round(time.monotonic()-t,2),'complete family members',g.shape[0]*g.shape[1]*g.shape[2],flush=True)

if __name__=='__main__': main()

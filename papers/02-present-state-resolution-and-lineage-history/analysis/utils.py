"""Package-owned utilities. Never write into the immutable previous tranche."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import h5py
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, svds
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent
PREVIOUS = Path(os.environ.get('RMMO_PREVIOUS', str(ROOT.parent/'previous_evidence'/'analysis')))
RESEARCH = PREVIOUS.parent
RAW = Path(os.environ.get('MELA_RAW_DIR', 'data/raw'))
CACHE = Path(os.environ.get('RMMO_REFEREE02_CACHE', 'data/referee02_cache'))
REPS = ['E7.5-R1', 'E7.5-R2', 'E7.5-R3']
PLAN = json.loads((ROOT/'preregistration/analysis_plan.json').read_text())
CONFIGS = PLAN['family']['configurations']
LAM = PLAN['family']['E1_lambda']
ALPHA = PLAN['family']['E1_alpha']
K = PLAN['family']['k']

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def stable(s): return int.from_bytes(hashlib.sha256(s.encode()).digest()[:8],'big')
def read(path):
    path=Path(path)
    with (gzip.open if path.suffix=='.gz' else open)(path,'rt',newline='') as f:
        return list(csv.DictReader(f))

def write(path, rows, fields=None):
    path=ROOT/path; path.parent.mkdir(parents=True,exist_ok=True)
    rows=list(rows); fields=fields or list(rows[0])
    with (gzip.open if path.suffix=='.gz' else open)(path,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n'); w.writeheader(); w.writerows(rows)

def js(path, value):
    p=ROOT/path; p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')

def unit(x):
    x=np.asarray(x,dtype=np.float64); n=np.linalg.norm(x,axis=1,keepdims=True)
    return x/np.where(n==0,1,n)

def scores(y,p):
    y,p=unit(y),unit(p)
    return {'NSE':float(np.mean(((y-p)**2).sum(1)/np.maximum((y*y).sum(1),1e-12))),
            'cosine':float(np.mean((y*p).sum(1))), 'MAE':float(np.mean(np.abs(y-p)))}

def distances(a,b):
    return np.maximum(0,(a*a).sum(1)[:,None]+(b*b).sum(1)[None,:]-2*a@b.T)

def standard(a,b):
    m=a.mean(0); s=a.std(0); s[s==0]=1
    return (a-m)/s,(b-m)/s

def history(h):
    h=h.copy(); h[:,:6]=np.log1p(h[:,:6]); return h

def load():
    old=RESEARCH/'relational_morphogenesis12'
    rows=read(old/'authority/transcriptome_lineage_join.csv.gz')
    hr={(r['replicate'],r['cellHash']):r for r in read(RESEARCH/'relational_morphogenesis11/devmap/lineage_history.csv.gz')}
    matrices=[]; hh=[]; genes=names=None
    for rep in REPS:
        rr=[r for r in rows if r['replicate']==rep]; ix=np.array([int(r['rowIndex']) for r in rr])
        with h5py.File(RAW/(rep+'.h5td'),'r') as h:
            x=sp.csr_matrix((h['X/data'][:],h['X/indices'][:],h['X/indptr'][:]),shape=tuple(h['X'].attrs['shape']))
            matrices.append(x[ix].astype(np.float64))
            dec=lambda a:np.array([v.decode() if isinstance(v,bytes) else str(v) for v in a])
            g,n=dec(h['var/gene_ids'][:]),dec(h['var/gene_name'][:])
            if genes is None: genes,names=g,n
            assert np.array_equal(g,genes) and np.array_equal(n,names)
            assert dec(h['obs/_index'][:])[ix].tolist()==[r['sourceCellId'] for r in rr]
            for r,i in zip(rr,ix):
                a=hr[(rep,r['cellHash'])]
                hh.append([float(a[k]) for k in ['depth','parentArity','grandparentArity','incomingEdits','cumulativeEdits','branchingAncestorCount']]+[float(h['obs/edit_frac'][i])])
    assert len(rows)==3704
    return rows,sp.vstack(matrices,format='csr'),genes,names,np.array(hh)

def panel_matrix(raw,ix):
    x=raw[:,ix].astype(np.float64); total=np.asarray(x.sum(1)).ravel()
    x=sp.diags(np.divide(10000.,total,out=np.zeros_like(total),where=total>0))@x
    x.data=np.log1p(x.data)
    return x.tocsr()

def project(x,genes,seed,dim):
    bucket=np.array([stable(f'{seed}|bucket|{g}')%dim for g in genes])
    sign=np.array([1. if stable(f'{seed}|sign|{g}')%2==0 else -1. for g in genes])
    return unit((x@sp.csr_matrix((sign,(np.arange(len(genes)),bucket)),shape=(len(genes),dim))).toarray())

def pca(matrix,tr,te):
    x=matrix[tr]; mu=np.asarray(x.mean(0)).ravel(); n,p=x.shape
    def mv(v): return np.asarray(x@np.asarray(v).reshape(-1))-mu@np.asarray(v).reshape(-1)
    def rmv(v): return np.asarray(x.T@np.asarray(v).reshape(-1))-mu*np.sum(v)
    def mm(v): return np.asarray(x@v)-mu@v
    def rmm(v): return np.asarray(x.T@v)-mu[:,None]*np.sum(v,axis=0)
    op=LinearOperator(x.shape,matvec=mv,rmatvec=rmv,matmat=mm,rmatmat=rmm,dtype=np.float64)
    rng=np.random.default_rng(stable('RMMO_P2_REFEREE02_TARGET_PCA_V0'))
    method='centered sparse ARPACK'
    try:
        _,s,vt=svds(op,k=32,tol=1e-8,maxiter=5000,v0=rng.normal(size=min(n,p)))
        order=np.argsort(s)[::-1]; s=s[order]; v=vt[order].T
    except Exception as error:
        method='training Gram fallback: '+type(error).__name__
        gram=(x@x.T).toarray(); xm=np.asarray(x@mu).ravel()
        gram-=xm[:,None]+xm[None,:]; gram+=float(mu@mu)
        val,u=eigh(gram,subset_by_index=(n-32,n-1)); order=np.argsort(val)[::-1]
        s=np.sqrt(np.maximum(val[order],1e-24)); u=u[:,order]
        v=(np.asarray(x.T@u)-mu[:,None]*u.sum(0))/s
    for j in range(32):
        if v[np.argmax(np.abs(v[:,j])),j]<0: v[:,j]*=-1
    a=np.asarray(x@v)-mu@v; b=np.asarray(matrix[te]@v)-mu@v
    return a,b,mu,v,s,method

class Bank:
    def __init__(self):
        self.rows,self.raw,self.genes,self.names,self.h=load()
        self.rep=np.array([r['replicate'] for r in self.rows]); self.s0=np.array([r['state'] for r in self.rows])
        self.cache={}; self.frames={}
    def matrices(self,reps):
        key=tuple(sorted(reps))
        if key in self.cache: return self.cache[key]
        tr=np.flatnonzero(np.isin(self.rep,key))
        det=np.array([np.asarray((self.raw[self.rep==r]>0).mean(0)).ravel() for r in key])
        eligible=np.all(det>=.01,axis=0)&(np.asarray(self.raw[tr].mean(0)).ravel()>0)
        eligible&=np.array([not str(n).startswith(('mt-','Mt-','MT-','Rpl','Rps')) for n in self.names])
        ix=np.flatnonzero(eligible); isS=np.array([stable('RMMO12_STATE_REALIZATION_V1|'+str(self.genes[j]))%10<7 for j in ix])
        s,t=ix[isS],ix[~isS]
        # Independently match the immutable published vocabulary, not just counts.
        old=read(PREVIOUS/'inductive_prediction/fold_panel_partitions.csv.gz')
        expected=[r for r in old if r['trainingEmbryos']=='+'.join(key)]
        assert s.tolist()==[int(r['geneIndex']) for r in expected if r['panel']=='state']
        assert t.tolist()==[int(r['geneIndex']) for r in expected if r['panel']=='target']
        self.cache[key]=(panel_matrix(self.raw,s),panel_matrix(self.raw,t),s,t)
        return self.cache[key]
    def representations(self,reps):
        key=tuple(sorted(reps))
        if key in self.frames: return self.frames[key]
        sm,tm,s,t=self.matrices(key); states={}; targets={}
        for name in CONFIGS:
            if name.startswith('target_seed'): continue
            if name=='pca32':
                old=np.load(PREVIOUS/('inductive_prediction/fitted_frames/pca_'+'_'.join(key)+'.npz'))
                assert np.array_equal(old['geneIndex'],s)
                states[name]=unit(np.asarray(sm@old['loadings'])-old['center']@old['loadings'])
            else:
                dim=int(name.split('_')[0][4:]); seed=int(name.split('seed')[1])
                label='RMMO12_STATE_HASH_PROJECTION_V1' if seed==0 else f'RMMO_P2_REFEREE01_STATE_SEED_{seed}'
                states[name]=project(sm,self.genes[s],label,dim)
        targets['canonical']=project(tm,self.genes[t],'RMMO12_REALIZATION_HASH_PROJECTION_V1',64)
        for i in range(1,8): targets[f'target_seed{i}']=project(tm,self.genes[t],f'RMMO_P2_REFEREE01_TARGET_SEED_{i}',64)
        self.frames[key]=(states,targets); return states,targets

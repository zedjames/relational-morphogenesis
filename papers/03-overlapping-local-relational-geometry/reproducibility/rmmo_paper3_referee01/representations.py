#!/usr/bin/env python3
"""Conventional outcome-blind PCA controls using the frozen feature partition."""
import datetime
import argparse
import os
import h5py
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, svds
from common import ROOT,REPS,check_freeze,load_state,parent,read,seed,sha,unit,js


def build_matrix(panel,rows):
    partition=[r for r in read(parent(12,'genes/state_realization_partition.csv.gz')) if r['panel']==panel]
    columns=np.asarray([int(r['geneIndex']) for r in partition])
    genes=np.asarray([r['geneId'] for r in partition])
    parts=[];raw_receipts=[]
    rawdir=os.environ.get('MELA_RAW_DIR','data/raw')
    for rep in REPS:
        source=__import__('pathlib').Path(rawdir)/(rep+'.h5td')
        selected=[r for r in rows if r['replicate']==rep]
        ix=np.asarray([int(r['rowIndex']) for r in selected])
        with h5py.File(source,'r') as f:
            decode=lambda a:np.asarray([v.decode() if isinstance(v,bytes) else str(v) for v in a])
            assert decode(f['obs/_index'][:])[ix].tolist()==[r['sourceCellId'] for r in selected]
            assert np.array_equal(decode(f['var/gene_ids'][:])[columns],genes)
            raw=sp.csr_matrix((f['X/data'][:],f['X/indices'][:],f['X/indptr'][:]),shape=tuple(f['X'].attrs['shape']))[ix].astype(np.float64)
            total=np.asarray(raw.sum(1)).ravel()
            x=sp.diags(np.divide(10000.,total,out=np.zeros_like(total),where=total>0))@raw[:,columns]
            x=x.tocsr();x.data=np.log1p(x.data);parts.append(x)
        raw_receipts.append({'embryo':rep,'filename':source.name,'bytes':source.stat().st_size,'sha256':sha(source)})
    matrix=sp.vstack(parts,format='csr')
    destination=ROOT/'data_authority'/(panel+'_expression.npz');destination.parent.mkdir(parents=True,exist_ok=True)
    sp.save_npz(destination,matrix)
    js('data_authority/'+panel+'_expression_receipt.json',{'shape':list(matrix.shape),'nonzero_count':matrix.nnz,'sha256':sha(destination),'raw_sources':raw_receipts,'normalization':'log1p(10000 * counts / whole-cell library)','feature_partition_sha256':sha(parent(12,'genes/state_realization_partition.csv.gz')),'annotation_history_fields_read':[],'panel':panel})
    return matrix


def pca(x,k,label):
    mu=np.asarray(x.mean(0)).ravel()
    def mv(v):v=np.asarray(v).reshape(-1);return np.asarray(x@v)-mu@v
    def rmv(v):v=np.asarray(v).reshape(-1);return np.asarray(x.T@v)-mu*v.sum()
    def mm(v):return np.asarray(x@v)-mu@v
    def rmm(v):return np.asarray(x.T@v)-mu[:,None]*np.sum(v,axis=0)
    op=LinearOperator(x.shape,matvec=mv,rmatvec=rmv,matmat=mm,rmatmat=rmm,dtype=np.float64)
    random=np.random.default_rng(seed('PCA|'+label))
    with np.errstate(all='ignore'):
        u,s,vt=svds(op,k=k,tol=1e-8,maxiter=20000,v0=random.normal(size=min(x.shape)))
        order=np.argsort(s)[::-1];s=s[order];v=vt[order].T
        for j in range(k):
            if v[np.argmax(np.abs(v[:,j])),j]<0:v[:,j]*=-1
        scores=np.asarray(x@v)-mu@v
    assert np.isfinite(scores).all() and np.isfinite(v).all()
    # Independent sparse row arithmetic cross-check, avoiding matrix BLAS path.
    manual=np.asarray([np.sum(x[i].data[:,None]*v[x[i].indices],axis=0)-np.sum(mu[:,None]*v,axis=0) for i in range(8)])
    assert np.max(np.abs(manual-scores[:8]))<1e-8
    residual=np.linalg.norm(op.matmat(v)-scores)/max(1.,np.linalg.norm(scores))
    assert residual<1e-8
    return scores,v,mu,s,float(residual)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--from-frozen-expression',action='store_true');args=ap.parse_args()
    check_freeze();rows,_=load_state()
    x=sp.load_npz(ROOT/'data_authority/state_expression.npz') if args.from_frozen_expression else build_matrix('state',rows)
    for dim in (32,64):
        label='pca'+str(dim);scores,v,mu,s,res=pca(x,dim,label)
        path=ROOT/'representation_families'/label/'representation.npz';path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path,embedding=unit(scores),scores=scores,loadings=v,mean=mu,singular_values=s)
        js('representation_families/'+label+'/representation_spec.json',{'dimensions':dim,'sha256':sha(path),'expression_sha256':sha(ROOT/'data_authority/state_expression.npz'),'centering':'all frozen source cells; no outcome/labels/history','solver':'centered sparse ARPACK','tol':1e-8,'maxiter':20000,'sign':'largest absolute loading positive','row_unit_normalized':True,'independent_row_arithmetic_max_tolerance':1e-8,'projection_residual':res,'target_panel_accessed_for_fit':False,'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
        print(label,'built and frozen',flush=True)
    target=sp.load_npz(ROOT/'data_authority/realization_expression.npz') if args.from_frozen_expression else build_matrix('realization',rows)
    scores,v,mu,s,res=pca(target,32,'target_pca32')
    path=ROOT/'representation_families/target_pca_control/representation.npz';path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,embedding=unit(scores),scores=scores,loadings=v,mean=mu,singular_values=s)
    js('representation_families/target_pca_control/representation_spec.json',{'dimensions':32,'sha256':sha(path),'expression_sha256':sha(ROOT/'data_authority/realization_expression.npz'),'state_geometry_changed':False,'fit_after_state_PCA_freeze':True,'projection_residual':res,'probe_authority':'held-out contemporaneous molecular panel, not independent replication'})
    print('target PCA32 built',flush=True)


if __name__=='__main__':main()

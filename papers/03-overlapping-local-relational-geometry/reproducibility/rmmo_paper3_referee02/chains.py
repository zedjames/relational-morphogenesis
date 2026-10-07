"""True proposal sweeps; persistent symmetric switch chains, not fresh starts."""
import hashlib
import ctypes
import subprocess
import tempfile
import os
from collections import defaultdict
from pathlib import Path
import numpy as np
from engine import ROOT,seed,sha
from rewire import make_rewirer

class Chain:
    def __init__(self,left,right,rows,bins,ixleft,ixright,group_by_source=True):
        self.obj=make_rewirer(left,right,rows,bins,ixleft,ixright)
        # Historical landmark kernels retain source-S0-conditioned target
        # degrees. R4 asks only for global target degrees and per-row target
        # S0/anchor-bin counts: its proposals must not retain that extra margin.
        if not group_by_source:
            groups=defaultdict(list)
            for i,(target,decile) in enumerate(zip(right,self.obj.edge_bins)):
                groups[(rows[int(target)]['state'],int(decile))].append(i)
            self.obj.groups=[np.asarray(g,dtype=np.int32) for g in groups.values() if len(g)>=2]
            self.obj.groupsizes=[len(g) for g in groups.values()]
        self.initial=self.obj.right.copy();self.current=self.initial.copy()
        self.proposals=0;self.accepted=0
        self.occupied=np.zeros(self.obj.bins.shape,dtype=np.uint8)
        self.occupied[self.obj.left,self.current]=1
        cache=Path(tempfile.gettempdir())/'rmmo_p3_referee02_kernel';cache.mkdir(exist_ok=True)
        binary=cache/('switch_'+sha(ROOT/'switch.cpp')[:16]+'.so')
        if not binary.exists():
            descriptor,temp_path=tempfile.mkstemp(prefix='switch_build_',suffix='.so',dir=cache);os.close(descriptor)
            subprocess.run(['c++','-std=c++17','-O3','-shared','-fPIC',str(ROOT/'switch.cpp'),'-o',temp_path],check=True)
            os.replace(temp_path,binary)
        library=ctypes.CDLL(str(binary));self.kernel=library.switch_block
        ptr=ctypes.c_void_p
        self.kernel.argtypes=[ptr,ptr,ptr,ptr,ctypes.c_int64,ptr,ctypes.c_int32,ptr,ptr]
        self.kernel.restype=ctypes.c_int64
        self.accepted_kernel=library.accepted_block
        self.accepted_kernel.argtypes=self.kernel.argtypes+[ctypes.c_int64,ptr]
        self.accepted_kernel.restype=ctypes.c_int64
        self.pending=None
        groups=[g for g in self.obj.groups];assigned={int(i) for g in groups for i in g}
        groups += [np.asarray([i],dtype=np.int32) for i in range(len(self.current)) if i not in assigned]
        self.flat=np.concatenate(groups).astype(np.int32) if groups else np.array([],np.int32)
        self.starts=np.zeros(len(self.current),dtype=np.int64);self.positions=self.starts.copy();self.sizes=self.starts.copy();offset=0
        for g in groups:
            self.starts[g]=offset;self.positions[g]=np.arange(len(g));self.sizes[g]=len(g);offset+=len(g)
    def reset(self):
        self.current=self.initial.copy();self.proposals=0;self.accepted=0
        self.occupied.fill(0);self.occupied[self.obj.left,self.current]=1
        self.pending=None
    def accepted_step(self,random,requested,max_attempts):
        accepted=attempted=0;obj=self.obj
        if len(self.current)==0:return {'accepted_switches':0,'attempted_switches':0,'complete':False}
        while accepted<requested and attempted<max_attempts:
            if self.pending is None:
                first=random.integers(len(self.current),size=min(8192,max_attempts-attempted))
                size=self.sizes[first];off=random.integers(np.maximum(size-1,1))
                off=np.where(size>1,off+(off>=self.positions[first]),0)
                second=self.flat[self.starts[first]+off]
                self.pending=(np.ascontiguousarray(first,dtype=np.int32),np.ascontiguousarray(second,dtype=np.int32))
            aa,bb=self.pending;used=np.zeros(1,dtype=np.int64)
            available=min(len(aa),max_attempts-attempted)
            gained=int(self.accepted_kernel(obj.left.ctypes.data,self.current.ctypes.data,aa.ctypes.data,bb.ctypes.data,
                available,obj.bins.ctypes.data,obj.bins.shape[1],obj.edge_bins.ctypes.data,self.occupied.ctypes.data,
                requested-accepted,used.ctypes.data))
            count=int(used[0]);accepted+=gained;attempted+=count
            self.pending=(aa[count:],bb[count:]) if count<len(aa) else None
        self.accepted+=accepted;self.proposals+=attempted;self.validate()
        return {'accepted_switches':accepted,'attempted_switches':attempted,'complete':accepted==requested}
    def block(self,first,second):
        obj=self.obj;aa=np.ascontiguousarray(first,dtype=np.int32);bb=np.ascontiguousarray(second,dtype=np.int32)
        return int(self.kernel(obj.left.ctypes.data,self.current.ctypes.data,aa.ctypes.data,bb.ctypes.data,
            len(aa),obj.bins.ctypes.data,obj.bins.shape[1],obj.edge_bins.ctypes.data,self.occupied.ctypes.data))
    def random_group_block(self,random,group,n):
        # Draw the entire first-index stream, then stream the second indices.
        # This keeps the original RNG order and exact switch order, without
        # simultaneously allocating four enormous proposal arrays at tier3.
        first=random.integers(len(group),size=n);accepted=0
        for start in range(0,n,65536):
            aa=first[start:start+65536]
            second=random.integers(len(group)-1,size=len(aa));second+=second>=aa
            accepted+=self.block(group[aa],group[second])
        return accepted
    def legacy_draw(self,random,sweeps=10):
        self.reset();attempted=0;accepted=0
        for group in self.obj.groups:
            n=min(max(2*len(group),20),5000)*sweeps
            accepted+=self.random_group_block(random,group,n);attempted+=n
        self.proposals=attempted;self.accepted=accepted;self.validate()
        return self.obj.right_nodes[self.current].copy(),{'accepted_swaps':accepted,'proposals':attempted}
    def step(self,random,sweeps=1):
        obj=self.obj;budget=int(sweeps*len(self.current))
        # Select first edge uniformly. Ineligible singleton strata are no-op
        # proposals. Within a stratum the second edge is uniform and distinct.
        sizes=np.asarray([len(g) for g in obj.groups],float)
        weights=np.r_[sizes,max(0,len(self.current)-sizes.sum())]/max(1,len(self.current))
        counts=random.multinomial(budget,weights)
        accepted=0
        for g,n in zip(obj.groups,counts):
            accepted+=self.random_group_block(random,g,int(n))
        self.proposals+=budget;self.accepted+=int(accepted)
        self.validate()
        return obj.right_nodes[self.current].copy(),{'attempted_switches':budget,'accepted_switches':int(accepted),
            'acceptance_fraction':float(accepted/budget) if budget else 0.,'edge_count':len(self.current),
            'graph_hash':self.graph_hash(),'edge_symmetric_difference_observed':self.difference(self.initial),
            'degree_mismatch':0,'S0_mismatch':0,'distance_bin_mismatch':0}
    def validate(self):
        obj=self.obj
        assert np.array_equal(np.bincount(self.current,minlength=len(obj.right_nodes)),obj.degree)
        assert np.array_equal(obj.bins[obj.left,self.current],obj.edge_bins)
        assert len(np.unique(obj.left.astype(np.int64)*len(obj.right_nodes)+self.current))==len(self.current)
        assert np.array_equal(obj.target_s0(obj.right_nodes[self.current]),obj.target_s0(obj.right_nodes[self.initial]))
    def difference(self,other):
        obj=self.obj;width=len(obj.right_nodes)
        a=obj.left.astype(np.int64)*width+self.current
        b=obj.left.astype(np.int64)*width+other
        return len(np.setxor1d(a,b))
    def graph_hash(self):
        obj=self.obj;width=len(obj.right_nodes)
        codes=np.sort(obj.left.astype(np.int64)*width+self.current)
        return hashlib.sha256(codes.tobytes()).hexdigest()

def autocorrelation(values):
    x=np.asarray(values,float);x=x-x.mean();den=x@x
    if den<=1e-20:return {'lag10':0.,'tau':None,'ESS':None,'constant_statistic':True}
    ac=[float(x[:-k]@x[k:]/den) for k in range(1,min(len(x)//2,64))]
    positive=[]
    for a in ac:
        if a<=0:break
        positive.append(a)
    tau=max(1.,1+2*sum(positive))
    return {'lag10':ac[9] if len(ac)>9 else None,'tau':tau,'ESS':float(len(x)/tau),'constant_statistic':False}

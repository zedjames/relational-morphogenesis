"""Symmetric stratified finite switch chains, with realized invariants checked."""
import ctypes
import subprocess
from collections import defaultdict
import numpy as np
from common import CACHE, ROOT, sha


class Rewirer:
    def __init__(self, source, target, rows, zbin, left_nodes, right_nodes):
        self.left_nodes, self.right_nodes = np.asarray(left_nodes), np.asarray(right_nodes)
        self.left = np.ascontiguousarray(np.searchsorted(left_nodes, source), dtype=np.int32)
        self.right = np.ascontiguousarray(np.searchsorted(right_nodes, target), dtype=np.int32)
        assert np.array_equal(left_nodes[self.left], source) and np.array_equal(right_nodes[self.right], target)
        self.bins = np.ascontiguousarray(zbin, dtype=np.uint8)
        self.edge_bins = self.bins[self.left, self.right].copy()
        groups = defaultdict(list)
        for i, (s, t, decile) in enumerate(zip(source, target, self.edge_bins)):
            groups[(rows[int(s)]['state'], rows[int(t)]['state'], int(decile))].append(i)
        self.groups = [np.asarray(pos, dtype=np.int32) for pos in groups.values() if len(pos) >= 2]
        self.groupsizes = [len(p) for p in groups.values()]
        CACHE.mkdir(parents=True, exist_ok=True)
        binary = CACHE / ('rewire_' + sha(ROOT / 'rewire.cpp')[:16] + '.so')
        if not binary.exists(): subprocess.run(['c++', '-std=c++17', '-O3', '-shared', '-fPIC', str(ROOT/'rewire.cpp'), '-o', str(binary)], check=True)
        self.lib = ctypes.CDLL(str(binary)).switches
        ptr = ctypes.c_void_p
        self.lib.argtypes = [ptr,ptr,ctypes.c_int64,ptr,ptr,ctypes.c_int64,ptr,ctypes.c_int32,ptr]
        self.lib.restype = ctypes.c_int64
        self.degree = np.bincount(self.right, minlength=len(right_nodes))

    def draw(self, random, sweeps=1):
        right = self.right.copy()
        first, second = [], []
        for pos in self.groups:
            budget = min(max(2*len(pos),20),5000)*sweeps
            a = random.integers(len(pos), size=budget)
            b = random.integers(len(pos)-1, size=budget)
            b += b >= a
            first.append(pos[a]); second.append(pos[b])
        aa = np.ascontiguousarray(np.concatenate(first) if first else [], dtype=np.int32)
        bb = np.ascontiguousarray(np.concatenate(second) if second else [], dtype=np.int32)
        accepted = self.lib(self.left.ctypes.data, right.ctypes.data, len(right), aa.ctypes.data, bb.ctypes.data, len(aa), self.bins.ctypes.data, self.bins.shape[1], self.edge_bins.ctypes.data)
        assert np.array_equal(np.bincount(right,minlength=len(self.right_nodes)), self.degree)
        assert np.array_equal(self.bins[self.left,right], self.edge_bins)
        assert len(np.unique(self.left.astype(np.int64)*len(self.right_nodes)+right)) == len(right)
        target = self.right_nodes[right]
        original = self.right_nodes[self.right]
        # Each proposal is within a source-S0/target-S0/bin group.
        assert np.array_equal(self.target_s0(target),self.target_s0(original))
        return target, {'accepted_swaps':int(accepted), 'proposals':len(aa), 'endpoint_displacement':float(np.mean(right!=self.right)) if len(right) else 0., 'degree_mismatch':0, 'distance_bin_mismatch':0, 'S0_mismatch':0}

    def target_s0(self, targets):
        # Bound after construction to avoid repeatedly copying the full carrier.
        return self.state[targets]


def make_rewirer(source, target, rows, bins, left_nodes, right_nodes):
    obj = Rewirer(source,target,rows,bins,left_nodes,right_nodes)
    obj.state = np.asarray([r['state'] for r in rows])
    return obj


def degree_only(source, target, random):
    """Uniform stub assignment conditional on simple graph; no fixed rejection cap."""
    trials = 0
    while True:
        new = random.permutation(target)
        trials += 1
        if len(set(zip(source.tolist(), new.tolist()))) == len(source):
            return new, {'stub_trials':trials, 'endpoint_displacement':float(np.mean(new!=target)), 'degree_mismatch':0}

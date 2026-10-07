"""Adversarial separation, sampler and tail controls (not biological evidence)."""
import itertools,unittest
import scipy.sparse as sp
from common05 import *
import construction05 as construct
def normalized(x):
    lib=np.asarray(x.sum(1)).ravel();x=(sp.diags(np.divide(10000.,lib,out=np.zeros_like(lib,dtype=float),where=lib>0))@x.astype(float)).tocsr();x.data=np.log1p(x.data);return x,lib
class Controls(unittest.TestCase):
    def test_zero_library(self):
        x,l=normalized(sp.csr_matrix([[0,0],[1,3]]));np.testing.assert_array_equal(l,[0,4]);self.assertEqual(x[0].nnz,0);self.assertTrue(np.isfinite(x.data).all())
    def test_opposite_panel_exclusion(self):
        x=np.array([[1,2,3,4],[0,0,7,8]]);y=x.copy();y[:,2:]*=99;np.testing.assert_array_equal(normalized(sp.csr_matrix(x[:,:2]))[0].toarray(),normalized(sp.csr_matrix(y[:,:2]))[0].toarray());self.assertFalse(np.array_equal(normalized(sp.csr_matrix(x))[0].toarray()[:,:2],normalized(sp.csr_matrix(y))[0].toarray()[:,:2]))
    def test_formula(self):
        np.testing.assert_allclose(normalized(sp.csr_matrix([[1,3],[0,0]]))[0].toarray(),[[np.log1p(2500),np.log1p(7500)],[0,0]],rtol=0,atol=1e-12)
    def test_forbidden_count_and_annotation_open(self):
        for rel in ('evaluation_inputs/heldout_expression.npz','authority/cell_carrier.csv.gz','authority/primary_GR.csv.gz'):
            with self.assertRaises(AssertionError):
                with construct.access_guard():open(ROOT/rel,'rb')
    def test_annotated_identity_rejected(self):
        with self.assertRaises(AssertionError):construct.compute(sp.csr_matrix((3704,14660)),[{'geneId':'x','geneIndex':'0'}],[{'cellHash':'x','replicate':'r','state':'forbidden'}],{},[])
    def test_tail_floor(self):
        d=family([3,3,3],np.ones((B,3)),False);self.assertEqual(d['sum_tail'],1/1025);self.assertEqual(d['min_tail'],1/1025)
    def test_sign_is_not_support(self):
        d=family([.9,.9,.9],np.linspace(0,2,B)[:,None]*np.ones((1,3)));self.assertTrue(all(r['relative_effect']>0 for r in d['components']));self.assertFalse(d['componentwise_family_supported'])
    def test_Tmin_is_not_marginal_maxT(self):
        d=json.loads((ROOT/'source_fidelity/q25/decision.json').read_text())['heldout_source_error'];self.assertLessEqual(d['min_tail'],.05);self.assertGreater(d['components'][0]['maxT_tail'],.05)
    def test_uniform_labeled_stub_pushforward(self):
        counts=Counter();left=[0,0,1,2]
        for p in itertools.permutations([0,0,1,2]):
            if len(set(zip(left,p)))==4:counts[tuple(sorted(zip(left,p)))]+=1
        self.assertGreater(len(counts),1);self.assertEqual(len(set(counts.values())),1)
    def test_subset_sampling_and_complement(self):
        for k in (2,3):
            a=Sampler.sample(np.arange(4),k,24000,np.random.default_rng(37));self.assertTrue(all(len(set(r))==k for r in a));c=Counter(map(tuple,np.sort(a,axis=1)));self.assertEqual(len(c),len(list(itertools.combinations(range(4),k))));self.assertLess(max(c.values())-min(c.values()),700)
    def test_frozen_counts(self):
        p=freeze();self.assertEqual(p['strict_subpanel_sizes'],[10300,10307,10268,10207]);self.assertEqual(p['consensus'],7);self.assertEqual(p['B'],1024);self.assertTrue(p['no_rescue_tranche'])
    def test_S1_fixed_before_sampling(self):
        a=json.loads((ROOT/'spatial_overlap/S1_construction_audit.json').read_text());self.assertTrue(a['fixed_before_sampling']);self.assertEqual(a['strata_sha256'],sha(ROOT/'spatial_overlap/fixed_cell_strata.npz'));self.assertFalse(a['S2_rerun'])
if __name__=='__main__':unittest.main(verbosity=2)

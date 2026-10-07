"""Post-reveal controls, object comparisons and transparent proposal accounting."""
import argparse,shutil
import scipy.sparse as sp
from common05 import *
def compare(a,b):
    a=set(map(int,a));b=set(map(int,b));return {'reference_count':len(a),'strict_count':len(b),'intersection':len(a&b),'gained':len(b-a),'lost':len(a-b),'Jaccard':len(a&b)/len(a|b) if a|b else 1.}
def prepare():
    revealed();path=ROOT/'geometry/hash_primary';path.mkdir(exist_ok=True);shutil.copyfile(P4/'geometry/hash_primary/pair_cores.npz',path/'pair_cores.npz')
    r,z,ix,c,g,o,m=arrays();p=np.load(P4/'geometry/hash/pair_cores.npz');oldc={n:(p[n+'_left'],p[n+'_right']) for n in PAIRS};og,oo=geometry(z,ix,oldc);records=[]
    for name in PAIRS:
        for label,j in [('source',0),('target',1)]:records.append({'object':label,'pair':name,**compare(oldc[name][j],c[name][j])})
        records.append({'object':'edge','pair':name,**compare(oldc[name][0]*3704+oldc[name][1],c[name][0]*3704+c[name][1])})
    for name,e in g.items():
        for q in SCALES:records.append({'object':'domain','orientation':name,'scale':q,**compare(og[name][q]['nodes'],e[q]['nodes'])})
    for name,_,_ in OVERLAPS:records.append({'object':'overlap','overlap':name,'scale':'q25',**compare(oo[name,'q25'],o[name,'q25'])})
    write('geometry/v4_vs_strict_object_identity.csv',records);print('Strict-v4 object comparison complete',flush=True)
def transparency():
    rr=read(ROOT/'source_fidelity/graph_audit.csv');out=[]
    for p in PAIRS:
        v=np.array([int(r['stub_trials']) for r in rr if r['pair']==p]);assert len(v)==B
        out.append({'pair':p,'attempted_stub_pairings':int(v.sum()),'simple_pairings':B,'acceptance_fraction':B/int(v.sum()),'rejected_nonsimple':int(v.sum()-B),'mean_proposals_per_accepted':float(v.mean()),'median_proposals_per_accepted':float(np.median(v)),'maximum_proposals':int(v.max()),'fresh_independent_pairing':True,'nonsimple_discarded_completely':True,'Markov_chain_used':False})
    write('source_fidelity/sampler_transparency.csv',out);(ROOT/'source_fidelity/degree_only_bank').mkdir(exist_ok=True)
    for rel in ('graph_audit.csv','degree_only_graphs.npz'):shutil.copyfile(ROOT/'source_fidelity'/rel,ROOT/'source_fidelity/degree_only_bank'/rel)
    md('source_fidelity/sampler_methods.md','Each null member uses a fresh independent uniformly permuted labeled target-stub pairing. Nonsimple proposals are discarded completely and resampled until a simple graph is obtained. No Markov chain or S0/distance restriction is used. Every simple graph with the fixed margins has the same multiplicity product(source_degree!)*product(target_degree!), so rejection induces the uniform simple-graph law. Acceptance fractions and all trial counts are reported per pair, not estimated from a finite walk. The same accepted pair draws are shared by the three scale sensitivities.');print('Degree sampler transparency complete',out,flush=True)
def retention_summary():
    for fn in ('component_statistics.csv','family_statistics.csv'):write('retention/'+fn,[{'level':f'R{i}',**r} for i in range(4) for r in read(ROOT/f'retention/R{i}/{fn}')])
    js('retention/decision.json',{f'R{i}':json.loads((ROOT/f'retention/R{i}/decision.json').read_text()) for i in range(4)})
def evaluation_raw():
    revealed();source=P4/'normalization/selected_eligible_counts.npz';raw=sp.load_npz(source);idx=np.load(P4/'normalization/count_index.npz');genes=read(ROOT/'authority/primary_GR.csv.gz');np.testing.assert_array_equal(idx['geneId'][idx['heldout_columns']],[g['geneId'] for g in genes]);sp.save_npz(ROOT/'evaluation_inputs/GR_raw_counts.npz',raw[:,idx['heldout_columns']].copy())
    js('evaluation_inputs/raw_GR_authority.json',{'created_utc':utc(),'copied_only_after_reveal':True,'source_selected_count_sha256':sha(source),'GR_raw_count_sha256':sha(ROOT/'evaluation_inputs/GR_raw_counts.npz'),'historical_hash_denominator_precision':'sum of raw count dtype before float64 multiplication; PCA expression sums float64 counts','geometry_changed':False})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','transparency','retention_summary','evaluation_raw']);a=p.parse_args();globals()[a.action]()

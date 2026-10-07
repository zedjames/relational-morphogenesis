"""One final normalization convention; immutable inferential vectors supply ranks."""
import argparse
from io08 import *
from calibration08 import catalog,inputs,S,P6,P7
def families():
    plan=json.loads((ROOT/'freeze/missing_calibration_plan.json').read_text());allowed={(s['bank'],s['metric']) for s in plan['banks']};specs=[f for f in catalog() if f['role']=='publication_primary' or (f['bank'],f['metric']) in allowed]
    assert len(specs)==36
    for f in specs:f['role']='publication_primary' if f['role']=='publication_primary' else 'manuscript_sensitivity'
    return specs
def calibration_source(f):
    old=S/'independent_calibration_bank'/f['bank']/'bank.npz'
    if old.exists():return old,f['offset'],'existing_Revision06'
    return ROOT/'independent_calibration_bank'/f['bank']/'bank.npz',0 if f['panel']=='state' else 3,'authorized_Final08'
def constants():
    assert not (ROOT/'statistical_authority/final_component_statistics.csv').exists(),'Fixed constants must precede final ranks';prior=read(S/'component_calibration.csv');rr=[];sources=[]
    for f in families():
        path,offset,owner=calibration_source(f);a=np.load(path);v=a['statistics'][:,offset:offset+3];_,rows,_,obs=inputs(f);np.testing.assert_array_equal(a['observed'][offset:offset+3],obs);assert v.shape==(1024,3) and np.isfinite(v).all();mu=v.mean(0);sd=v.std(0,ddof=1)
        pc=[r for r in prior if r['bank']==f['bank'] and r['panel']==f['panel'] and r['metric']==f['metric']]
        if pc:np.testing.assert_array_equal(mu,[float(x['mu_cal']) for x in pc]);np.testing.assert_array_equal(sd,[float(x['sigma_cal']) for x in pc])
        rr.extend({**f,'component':rows[i]['component'],'component_index':i,'mu_cal':float(mu[i]),'sigma_cal':float(sd[i]),'B_cal':1024,'ddof':1,'eta':1e-8,'calibration_source':owner,'calibration_offset':offset,'calibration_bank_sha256':sha(path)} for i in range(3));sources.append(dict(bank=f['bank'],panel=f['panel'],metric=f['metric'],source=owner,sha256=sha(path)))
    write('statistical_authority/component_calibration_final.csv',rr);js('freeze/component_constants_freeze.json',dict(constants_sha256=sha(ROOT/'statistical_authority/component_calibration_final.csv'),plan_sha256=sha(ROOT/'freeze/missing_calibration_plan.json'),utc=utc(),components=108,families=36,sources=sources,rescoring_not_yet_run=True,new_banks=7,existing_inferential_banks_unchanged=True));print('108 independent constants fixed; inferential ranks not inspected',flush=True)
def rank(v,o):return float((1+np.count_nonzero(v>=o))/1025)
def rescore():
    receipt=json.loads((ROOT/'freeze/component_constants_freeze.json').read_text());assert sha(ROOT/'statistical_authority/component_calibration_final.csv')==receipt['constants_sha256'];assert (ROOT/'freeze/constants_commit_receipt.json').exists(),'Commit constants before final ranks';cc=read(ROOT/'statistical_authority/component_calibration_final.csv');comp=[];fam=[];comparison=[];family_comparison=[]
    for f in families():
        path,rows,v,obs=inputs(f);s=f['direction'];c=[r for r in cc if all(str(f[k])==r[k] for k in ('bank','panel','metric'))];mu=np.array([float(r['mu_cal']) for r in c]);sd=np.array([float(r['sigma_cal']) for r in c]);z=s*(v-mu)/(sd+1e-8);zo=s*(obs-mu)/(sd+1e-8)
        original=next(r for r in read(path/'family_statistics.csv') if not r.get('panel') or r['panel']==f['panel'] and r['metric']==f['metric']);ff={**f,'T_sum':float(zo.sum()),'T_min':float(zo.min()),'sum_tail':rank(z.sum(1),zo.sum()),'min_tail':rank(z.min(1),zo.min()),'standardization':'independent_calibration'};fam.append(ff)
        for stat in ('sum','min'):
            tail0=float(original[stat+'_tail']);tail1=ff[stat+'_tail'];family_comparison.append({**f,'statistic':stat,'original_tail':tail0,'calibrated_tail':tail1,'original_supported':tail0<=.05,'calibrated_supported':tail1<=.05,'decision_changed':(tail0<=.05)!=(tail1<=.05)})
        for i,r in enumerate(rows):
            raw=rank(s*v[:,i],s*obs[i]);assert raw==float(r['conditional_tail']);tail=rank(z.max(1),zo[i]);old=float(r['maxT_tail']);comp.append({**f,**r,'panel':f['panel'],'metric':f['metric'],'original_Z':r['Z'],'Z':float(zo[i]),'maxT_tail':tail,'original_maxT_tail':old,'calibration_mu':float(mu[i]),'calibration_sigma':float(sd[i]),'standardization':'independent_calibration','T_sum':ff['T_sum'],'T_min':ff['T_min'],'sum_tail':ff['sum_tail'],'min_tail':ff['min_tail']});comparison.append({**f,'component':r['component'],'original_maxT_tail':old,'calibrated_maxT_tail':tail,'original_supported':old<=.05,'calibrated_supported':tail<=.05,'decision_changed':(old<=.05)!=(tail<=.05),'raw_tail_unchanged':True})
    assert len(comp)==len(comparison)==108 and len(fam)==36 and len(family_comparison)==72
    write('statistical_authority/final_component_statistics.csv',comp);write('statistical_authority/final_family_statistics.csv',fam);write('statistical_authority/final_maxT_comparison.csv',comparison);write('statistical_authority/final_family_decision_comparison.csv',family_comparison)
    changed=[r for r in comparison if r['decision_changed']];fc=[r for r in family_comparison if r['decision_changed']];js('statistical_authority/final_calibration_decision.json',dict(status='CALIBRATION_PRESERVED' if not changed and not fc else 'CLAIMS_REQUIRE_EXACT_UPDATE',families=36,components=108,component_decisions_changed=changed,family_events=72,family_decisions_changed=fc,all_reported_families_independently_calibrated=True,new_calibration_banks=7,old_inferential_banks_changed=0,target_PCA_reported=False,no_rescue=True,publication_freeze_complete=False))
    write('statistical_authority/final_scale_sensitivity.csv',[r for r in comp if r['panel']=='heldout' and r['scale'] in SCALES and r['bank'].split('/')[0] in ('matched_coherence','source_fidelity')]);write('statistical_authority/final_retention_effect_surface.csv',[dict(level=r['bank'].split('/')[-1],**r) for r in comp if r['bank'].startswith('retention/')]);print('Final calibrated ranks:',len(changed),'component decision changes;',len(fc),'family decision changes',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['constants','rescore']);globals()[p.parse_args().action]()

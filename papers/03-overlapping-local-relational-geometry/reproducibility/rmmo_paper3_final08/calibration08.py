"""Pure arithmetic on existing banks; never draws samples or reconstructs geometry."""
from io08 import *
P5=ROOT.parent/'rmmo_paper3_referee05';P6=ROOT.parent/'rmmo_paper3_revision06';P7=ROOT.parent/'rmmo_paper3_acceptance07';S=P6/'statistical_robustness'
def catalog():
    checks=json.loads((P7/'display_authority/maxT_numeric_receipt.json').read_text())['checks'];plan=json.loads((S/'calibration_plan.json').read_text());primary={(f['bank'],f['panel'],f['metric']) for f in plan['families']};primary.update((f'retention/R{i}',p,k) for i in range(4) for p in ('state','heldout') for k in ('defect','source_error'));out=[]
    for check in checks:
        b=check['bank'];owner='rmmo_paper3_acceptance07' if b.startswith('own/') else 'rmmo_paper3_referee05';rel=b[4:] if b.startswith('own/') else b.split('/',1)[1];path=ROOT.parent/owner/rel;rr=read(path/'component_statistics.csv');panel=check['panel'];metric=check['metric']
        if rel=='spatial_overlap/S1':panel,metric='spatial','overlap_size'
        chosen=[(i,r) for i,r in enumerate(rr) if not r.get('panel') or r['panel']==panel and r['metric']==metric];assert len(chosen)==3
        off=chosen[0][0];assert [i for i,_ in chosen]==list(range(off,off+3))
        scale=rel.rsplit('/',1)[-1] if rel.rsplit('/',1)[-1] in SCALES else 'q25' if rel.startswith(('retention/','representation/','spatial_overlap/','acceptance07/')) else 'q50_relation'
        is_primary=(rel,panel,metric) in primary;role='publication_primary' if is_primary else 'manuscript_sensitivity' if (rel.startswith(('representation/pca','acceptance07/')) and metric=='defect' or scale in ('q10','q50') and panel=='heldout' and (rel.startswith('matched_coherence/') and metric=='defect' or rel.startswith('source_fidelity/') and metric=='source_error')) else 'companion_provenance'
        cal=S/'independent_calibration_bank'/rel/'bank.npz';out.append(dict(owner=owner,bank=rel,family=rel,panel=panel,metric=metric,scale=scale,offset=off,direction=1 if rel=='spatial_overlap/S1' else -1,role=role,independent_bank_available=cal.exists()))
    assert len(out)==56 and len({(f['bank'],f['panel'],f['metric']) for f in out})==56;assert sum(f['role']=='publication_primary' for f in out)==22
    return out
def inputs(f):
    path=ROOT.parent/f['owner']/f['bank'];a=np.load(path/'bank.npz');off=f['offset'];rr=read(path/'component_statistics.csv')[off:off+3];v=a['statistics'][:,off:off+3];obs=a['observed'][off:off+3];assert v.shape==(1024,3) and np.isfinite(v).all() and np.isfinite(obs).all();return path,rr,v,obs
def rank(v,obs):return float((1+np.count_nonzero(v>=obs))/1025)
def main():
    assert not (ROOT/'freeze/missing_calibration_plan.json').exists(),'Historical pre-authorization diagnostic only; use publication_stats08.py and verify08.py for final authority'
    specs=catalog();prior=read(S/'component_calibration.csv');freeze=json.loads((S/'constants_freeze.json').read_text());assert sha(S/'component_calibration.csv')==freeze['constants_sha256'];comparisons=[];constants=[];family_results=[];input_records=[];prior_verified=0
    for f in specs:
        path,rows,v,obs=inputs(f);s=f['direction'];mu_old=v.mean(0);sd_old=v.std(0,ddof=1);z_old=s*(v-mu_old)/(sd_old+1e-8);zo_old=s*(obs-mu_old)/(sd_old+1e-8);mu=sd=z=zo=None
        for name in ('bank.npz','component_statistics.csv','family_statistics.csv'):input_records.append(dict(parent=f['owner'],artifact=f['bank']+'/'+name,sha256=sha(path/name)))
        if f['independent_bank_available']:
            calpath=S/'independent_calibration_bank'/f['bank']/'bank.npz';a=np.load(calpath);off=f['offset'];cv=a['statistics'][:,off:off+3];np.testing.assert_array_equal(a['observed'][off:off+3],obs);assert cv.shape==(1024,3);mu=cv.mean(0);sd=cv.std(0,ddof=1);z=s*(v-mu)/(sd+1e-8);zo=s*(obs-mu)/(sd+1e-8)
            pc=[r for r in prior if r['bank']==f['bank'] and r['panel']==f['panel'] and r['metric']==f['metric']]
            if pc:
                assert len(pc)==3;np.testing.assert_array_equal([float(r['mu_cal']) for r in pc],mu);np.testing.assert_array_equal([float(r['sigma_cal']) for r in pc],sd);prior_verified+=3
            for i in range(3):constants.append({**f,'component_index':i,'component':rows[i]['component'],'mu_cal':float(mu[i]),'sigma_cal':float(sd[i]),'constants_source':'existing_frozen_constants' if pc else 'already_scored_immutable_bank_columns','calibration_bank_sha256':sha(calpath)})
            family_results.append({**f,'T_sum':float(zo.sum()),'T_min':float(zo.min()),'sum_tail':rank(z.sum(1),zo.sum()),'min_tail':rank(z.min(1),zo.min()),'status':'COMPUTED_FROM_EXISTING_AUTHORITY'})
        else:family_results.append({**f,'status':'MISSING_MATCHING_INDEPENDENT_CALIBRATION'})
        for i,r in enumerate(rows):
            raw=rank(s*v[:,i],s*obs[i]);orig=rank(z_old.max(1),zo_old[i]);assert raw==float(r['conditional_tail']) and orig==float(r['maxT_tail']);new=rank(z.max(1),zo[i]) if z is not None else None
            comparisons.append({**f,'component':r['component'],'observed':float(obs[i]),'primitive_raw_tail':raw,'original_maxT_tail':orig,'calibrated_maxT_tail':new if new is not None else '', 'original_supported':orig<=.05,'calibrated_supported':new<=.05 if new is not None else 'NOT_COMPUTED','decision_changed':(orig<=.05)!=(new<=.05) if new is not None else 'NOT_COMPUTED','calibrated_Z':float(zo[i]) if zo is not None else '', 'status':'COMPUTED_FROM_EXISTING_AUTHORITY' if new is not None else 'MISSING_MATCHING_INDEPENDENT_CALIBRATION'})
    assert prior_verified==66;available=[f for f in specs if f['independent_bank_available']];missing=[f for f in specs if not f['independent_bank_available']];changed=[r for r in comparisons if r['decision_changed'] is True];primary=[r for r in comparisons if r['role']=='publication_primary'];assert len(primary)==66 and all(r['status']=='COMPUTED_FROM_EXISTING_AUTHORITY' for r in primary)
    write('statistical_authority/family_catalog.csv',specs);write('statistical_authority/component_calibration_available.csv',constants);write('statistical_authority/maxT_comparison.csv',comparisons);write('statistical_authority/family_statistics_available.csv',family_results);write('statistical_authority/input_hashes.csv',input_records)
    js('statistical_authority/coverage_decision.json',dict(status='INCOMPLETE_AWAITING_CALIBRATION_SCOPE_CHOICE',retained_families=len(specs),manuscript_families=sum(f['role']!='companion_provenance' for f in specs),primary_families=22,previous_frozen_components_verified=66,available_families=len(available),available_components=3*len(available),missing_families=len(missing),missing_components=3*len(missing),changed_available_components=changed,changed_primary_components=[r for r in primary if r['decision_changed']],no_new_random_bank=True,no_primitive_tails_changed=True,all_publication_component_decisions_compared=False,publication_freeze_complete=False))
    md('completion/calibration_coverage_issue.md','''# Final08 calibration authority issue — publication freeze incomplete

The prior 66 frozen constants cover 22 publication-primary families. All are
verified against the existing independent calibration draws; none of their
component support decisions changes. The complete retained-component catalog
contains 56 families (168 components), including companion-only provenance.

Four additional companion-only families can use already-scored columns of the
same immutable q25 calibration banks, without any new draw or geometry work.
Thus 26 retained families are calculable from existing authority; 30 are not.
This does not expand the primary manuscript claims or introduce a new test.

Missing corresponding banks/constants affect the q10/q50 matched and degree
sensitivities and hash/PCA controls. No other scale, representation or null law
has been substituted. Missing comparison cells are explicitly NOT_COMPUTED,
never reported as unchanged decisions. The 32 manuscript-facing families
include ten uncovered sensitivity families; 22 are covered.

Required choice remains: preserve the no-new-bank instruction and explicitly
retain historical standardization for uncovered sensitivities, or authorize
the missing matching calibration work. The full unified publication freeze
cannot be claimed under the original mutually incompatible requirements.

All original primitive tails and inferential banks remain unchanged. No new
random sampling, biological reconstruction, or null law has been performed.
Prepared display/manuscript changes remain provisional until this is resolved.
''')
    print('Existing-authority maxT arithmetic complete:',len(available),'families;',len(missing),'missing; changed primary components:',len([r for r in primary if r['decision_changed']]),flush=True)
if __name__=='__main__':main()

"""Strictly training-fitted non-hash target outcomes, with nested k selection."""
from utils import *

def main():
    bank=Bank(); frames={}; authorities=[]; tuning=[]; tables=[]; comparisons=[]; raw_tables=[]
    def frame(fit):
        key=tuple(sorted(fit))
        if key in frames: return frames[key]
        sm,tm,s,t=bank.matrices(key); tr=np.flatnonzero(np.isin(bank.rep,key)); other=np.flatnonzero(~np.isin(bank.rep,key))
        a,b,mu,v,singular,method=pca(tm,tr,other)
        raw=np.empty((3704,32)); raw[tr]=a; raw[other]=b; y=unit(raw)
        z=project(sm,bank.genes[s],'RMMO12_STATE_HASH_PROJECTION_V1',64)
        p=ROOT/'target_pca/fitted_frames'; p.mkdir(parents=True,exist_ok=True)
        filename='target_'+'_'.join(key)+'.npz'
        np.savez_compressed(p/filename,targetGeneIndex=t,stateGeneIndex=s,center=mu,loadings=v,singularValues=singular,trainIndices=tr,excludedIndices=other,rawCoordinates=raw,normalizedCoordinates=y)
        authorities.append({'fitEmbryos':list(key),'fitIndices':tr.tolist(),'excludedIndices':other.tolist(),'fitCellCount':len(tr),'targetGenes':len(t),'stateGenes':len(s),'targetDimension':32,'method':method,'frame':filename,'frameSHA256':digest(p/filename),'targetPanelOwnDenominator':True,'panelsDisjoint':not bool(set(s)&set(t)),'signRule':'largest absolute loading positive; first index breaks ties','rawCoordinateFitMeanMaxAbs':float(np.max(np.abs(a.mean(0))))})
        frames[key]=(z,y,raw,tr); print('Fitted target-only TRAIN frame',key,flush=True)
        return frames[key]
    for held in REPS:
        fit=[r for r in REPS if r!=held]; losses={k:[] for k in K}
        for validation in fit:
            inner=[r for r in fit if r!=validation]; z,y,raw,tr=frame(inner); te=np.flatnonzero(bank.rep==validation)
            neighbors=np.argsort(distances(z[te],z[tr]),axis=1,kind='stable')[:,:50]
            cumulative=np.cumsum(y[tr][neighbors],axis=1)
            for k in K: losses[k].append(scores(y[te],cumulative[:,k-1]/k)['NSE'])
        k=min(K,key=lambda k:(float(np.mean(losses[k])),k))
        tuning.append({'heldOutEmbryo':held,'outerTrainingEmbryos':fit,'grid':[{'k':a,'innerNSE':losses[a],'meanInnerNSE':float(np.mean(losses[a]))} for a in K],'selectedK':k})
        z,y,raw,tr=frame(fit); te=np.flatnonzero(bank.rep==held)
        nb=np.argsort(distances(z[te],z[tr]),axis=1,kind='stable')[:,:k]
        globalmean=y[tr].mean(0); pred0=np.tile(globalmean,(len(te),1)); raw0=np.zeros((len(te),32)); s0=pred0.copy(); rawS0=raw0.copy(); fallback=0
        for j,i in enumerate(te):
            sel=tr[bank.s0[tr]==bank.s0[i]]
            if len(sel): s0[j]=y[sel].mean(0); rawS0[j]=raw[sel].mean(0)
            else: fallback+=1
        pred={'molecular':y[tr][nb].mean(1),'S0Only':s0,'globalMean':pred0}
        pred_raw={'molecular':raw[tr][nb].mean(1),'S0Only':rawS0,'globalMean':raw0}
        met={model:scores(y[te],p) for model,p in pred.items()}
        variance=float(np.mean(np.sum(raw[tr]**2,axis=1)))
        for model in pred:
            tables.append({'held_out_embryo':held,'model':model,'selected_k':k,'test_cells':len(te),'fallback_cells':fallback if model=='S0Only' else 0,**met[model]})
            raw_tables.append({'held_out_embryo':held,'model':model,'training_PC_variance':variance,'training_scaled_PC_error':float(np.mean(np.sum((raw[te]-pred_raw[model])**2,axis=1))/variance)})
        gm=met['S0Only']['NSE']-met['molecular']['NSE']; gg=met['globalMean']['NSE']-met['molecular']['NSE']
        comparisons.append({'held_out_embryo':held,'selected_k':k,'molecular_NSE':met['molecular']['NSE'],'S0_NSE':met['S0Only']['NSE'],'global_mean_NSE':met['globalMean']['NSE'],'gain_over_S0':gm,'fractional_gain_over_S0':gm/met['S0Only']['NSE'],'gain_over_global':gg,'cosine':met['molecular']['cosine'],'MAE':met['molecular']['MAE']})
        d=ROOT/'target_pca/predictions'; d.mkdir(exist_ok=True)
        np.savez_compressed(d/(held+'.npz'),trainIndices=tr,testIndices=te,neighborTrainIndices=tr[nb],actual=y[te],actualRaw=raw[te],**pred,**{m+'Raw':p for m,p in pred_raw.items()})
        print('Target PCA outer result',comparisons[-1],flush=True)
    write('target_pca/foldwise_prediction.csv',tables)
    write('target_pca/baseline_comparison.csv',comparisons)
    write('target_pca/molecular_increment_over_s0.csv',comparisons)
    write('target_pca/unnormalized_PC_sensitivity.csv',raw_tables)
    js('target_pca/fold_target_authority.json',authorities); js('target_pca/nested_tuning_receipts.json',tuning)
    old=read(PREVIOUS/'inductive_prediction/prediction_hash64_seed0.csv'); matched=[]
    for row in comparisons:
        held=row['held_out_embryo']; a={r['model']:float(r['NSE']) for r in old if r['heldOutEmbryo']==held}
        matched += [{'target_family':'signed_hash64','held_out_embryo':held,'molecular_NSE':a['knn'],'S0_NSE':a['S0Only'],'matched_gain':a['S0Only']-a['knn'],'fractional_gain':(a['S0Only']-a['knn'])/a['S0Only']},{'target_family':'training_fitted_PCA32','held_out_embryo':held,'molecular_NSE':row['molecular_NSE'],'S0_NSE':row['S0_NSE'],'matched_gain':row['gain_over_S0'],'fractional_gain':row['fractional_gain_over_S0']}]
    write('target_pca/target_representation_comparison.csv',matched)
    positive=sum(r['gain_over_S0']>0 for r in comparisons)
    js('target_pca/decision.json',{'status':'NONHASH_TARGET_REPLICATES_3OF3' if positive==3 else 'NONHASH_TARGET_PARTIAL' if positive else 'NONHASH_TARGET_FAILS','positiveFolds':positive,'globalMeanPositiveFolds':sum(r['gain_over_global']>0 for r in comparisons),'foldwise':comparisons,'populationReplication':False,'representationIndependent':False,'targetDimension':32,'heldOutEmbryoExcludedFromAllFits':True})

if __name__=='__main__': main()

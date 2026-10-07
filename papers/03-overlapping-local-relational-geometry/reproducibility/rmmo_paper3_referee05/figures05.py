"""Six actual, source-CSV-backed vector review figures; no schematic placeholders."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common05 import *
from report05 import select,decision,TEST_NOTE
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'figure.dpi':120})
def frame(n,title):
    f,axes=plt.subplots(2,2,figsize=(11.2,8.0));f.suptitle(f'Figure {n}. {title}',fontsize=15);return f,axes.ravel(),[]
def finish(n,f,source):
    f.tight_layout(rect=(0,0,1,.94));p=ROOT/'figures';p.mkdir(exist_ok=True);f.savefig(p/f'figure{n:02}.pdf');f.savefig(p/f'figure{n:02}.png',dpi=160);plt.close(f);write(f'figures/figure{n:02}_authority.csv',source)
def bars(a,title,labels,vals,ylabel):a.bar(labels,vals,color='#397b99');a.set(title=title,ylabel=ylabel);a.tick_params(axis='x',labelrotation=25)
def text(a,title,s):a.axis('off');a.set_title(title,loc='left');a.text(0,.93,s,va='top',wrap=True,fontsize=10,transform=a.transAxes)
def main():
    f,a,s=frame(1,'Strict construction feature separation');lanes=read(ROOT/'strict_lanes/lane_registry.csv');gs=read(ROOT/'strict_lanes/lane_intersection_audit.csv');obj=read(ROOT/'geometry/v4_vs_strict_object_identity.csv')
    bars(a[0],'A  Frozen nine-lane bank',[r['lane'] for r in lanes],[int(r['genes']) for r in lanes],'Own-lane gene count');s.extend({'panel':'A',**r} for r in lanes)
    text(a[1],'B  Construction boundary','S0: full GS, primary projection\nS1–S4: full GS, frozen projections\nS5–S8: historical 70/30 assignment,\nrestricted to GS before construction.\nAll denominators: own selected genes only.\nOnly cell IDs and embryo IDs supplied.');s.extend({'panel':'B',**r} for r in read(ROOT/'freeze/strict_lane_plan.csv'))
    bars(a[2],'C  Complete union audit',['GS union','GR intersection'],[14660,0],'Gene identities');s.extend([{'panel':'C','quantity':'union_genes','value':14660},{'panel':'C','quantity':'union_intersection_GR','value':0}])
    rr=[r for r in obj if r['object']=='edge' or r['object']=='domain' and r['scale']=='q25'];bars(a[3],'D  Referee04 → strict finite objects',[r['pair'] or r['orientation'] for r in rr],[float(r['Jaccard']) for r in rr],'Jaccard');a[3].set_ylim(0,1);s.extend({'panel':'D',**r} for r in rr);finish(1,f,s)
    f,a,s=frame(2,'Construction-held-out molecular validation')
    specs=[('A  Landmark-edge GR agreement','edge_agreement',None,None),('B  Matched GR coherence','matched_coherence/q25','heldout','defect'),('C  Degree-controlled GR fidelity','source_fidelity/q25','heldout','source_error')]
    for ax,(title,rel,p,k) in zip(a,specs):
        rr=select(rel,p,k);x=np.arange(3);ax.bar(x-.16,[float(r['observed']) for r in rr],.32,label='Observed');ax.bar(x+.16,[float(r['null_mean']) for r in rr],.32,yerr=[float(r['null_sd']) for r in rr],capsize=3,label='Null mean ± draw SD');ax.set_xticks(x,[r['component'] for r in rr]);ax.set(title=title,ylabel='Molecular chord error/defect');ax.legend(fontsize=7);s.extend({'panel':title[0],**r} for r in rr)
    fi=decision('source_fidelity/q25','heldout','source_error');text(a[3],'D  Different family/marginal rejection events',f"All three fidelity reductions positive.\nSum tail = {fi['sum_tail']:.6g}\nMinimum tail = {fi['min_tail']:.6g}\nO1 maxT tail = {fi['components'][0]['maxT_tail']:.6g} (unsupported)\n\nT_min compares with the null MINIMUM.\nComponent maxT compares with null MAXIMUM.\nFamily support does not imply all marginal support.\nDraw SD is not biological-replicate uncertainty.");s.append({'panel':'D',**{k:v for k,v in fi.items() if k!='components'}});finish(2,f,s)
    f,a,s=frame(3,'Spatial and descriptive biological geometry');rr=select('spatial_overlap/S1');bars(a[0],'A  Strict q25 overlaps',[r['component'] for r in rr],[float(r['observed']) for r in rr],'Shared source cells');s.extend({'panel':'A',**r} for r in rr)
    bank=np.load(ROOT/'spatial_overlap/S1/bank.npz');a[1].boxplot(bank['statistics'],tick_labels=['O1','O2','O3'],showfliers=False);a[1].scatter(np.arange(1,4),bank['observed'],color='red',label='Observed',zorder=3);a[1].set(title='B  Fixed-stratum S1 reference (B=1024)',ylabel='Overlap cells');a[1].legend(fontsize=8);s.extend({'panel':'B','draw':b,'component':f'O{i+1}','overlap_cells':v} for b,row in enumerate(bank['statistics']) for i,v in enumerate(row))
    for ax,label,file,objects in [(a[2],'C  Landmark germ-layer composition','landmark_annotation_summary',['P12','P13','P23']),(a[3],'D  Overlap germ-layer composition','overlap_annotation_summary',['O1','O2','O3'])]:
        rr=[r for r in read(ROOT/'biology'/f'{file}.csv') if r['normalization']=='isolated' and r['field']=='germLayer' and r['object'] in objects];cats=sorted({r['category'] for r in rr});bottom=np.zeros(3)
        for cat in cats:
            v=[next(float(r['fraction']) for r in rr if r['object']==o and r['category']==cat) for o in objects];ax.bar(objects,v,bottom=bottom,label=cat);bottom+=v
        ax.set(title=label,ylabel='Descriptive fraction');ax.legend(fontsize=6,ncol=2);s.extend({'panel':label[0],**r} for r in rr)
    finish(3,f,s)
    f,a,s=frame(4,'Competitive baselines on declared common support');base=read(ROOT/'baselines/comparison.csv')
    for ax,label,b in zip(a[:3],('A  Matched-size kNN','B  Native MNN5','C  Eligible fixed-size MNN5'),('matched_knn','mnn_k5','mnn_fixed_size')):
        rr=[r for r in base if r['baseline']==b and r['scale']=='q25' and r['orientation'].startswith('O')];x=np.arange(3)
        for j,p in enumerate(('state','heldout')):ax.bar(x+(j-.5)*.32,[float(next(r for r in rr if r['panel']==p and r['orientation']==f'O{i+1}')['paired_difference_mean']) for i in range(3)],.32,label=p)
        ax.axhline(0,color='black',lw=.6);ax.set_xticks(x,['O1','O2','O3']);ax.set(title=label,ylabel='Native − baseline error (negative favors native)');ax.legend(fontsize=8);s.extend({'panel':label[0],**r} for r in rr)
    support=read(ROOT/'baselines/support.csv');rr=[r for r in base if r['baseline']=='mnn_fixed_size' and r['scale']=='q25' and r['panel']=='heldout' and r['orientation'].startswith('O')];bars(a[3],'D  Fixed-size MNN common overlap support',[r['orientation'] for r in rr],[int(r['cells']) for r in rr],'Eligible/common source cells');s.extend({'panel':'D',**r} for r in support);s.extend({'panel':'D_summary',**r} for r in rr);finish(4,f,s)
    f,a,s=frame(5,'Tested representation persistence, not object equality')
    for ax,label,v in zip(a[:3],('A  Strict hash','B  PCA32 state control','C  PCA64 state control'),('hash','pca32','pca64')):
        rel='matched_coherence/q25' if v=='hash' else 'representation/'+v;rr=select(rel,'heldout','defect');bars(ax,label,[r['component'] for r in rr],[100*float(r['relative_effect']) for r in rr],'Matched GR defect reduction (%)');s.extend({'panel':label[0],**r} for r in rr)
    rr=[r for r in read(ROOT/'representation/object_identity.csv') if r['reference']=='hash'];kinds=['source','target','edge','domain','overlap']
    for j,v in enumerate(('pca32','pca64')):
        vals=[np.mean([float(r['Jaccard']) for r in rr if r['object']==k and r['comparison']==v]) for k in kinds];a[3].bar(np.arange(5)+(j-.5)*.32,vals,.32,label=v)
    a[3].set_xticks(np.arange(5),kinds);a[3].set(title='D  Strict consensus vs single-lane PCA',ylabel='Mean object Jaccard',ylim=(0,1));a[3].legend(fontsize=8);s.extend({'panel':'D',**r} for r in read(ROOT/'representation/object_identity.csv'));finish(5,f,s)
    f,a,s=frame(6,'Final within-carrier authority');rr=read(ROOT/'retention/family_statistics.csv');v=np.array([[float(next(r for r in rr if r['level']==f'R{i}' and r['panel']==p and r['metric']==m)['min_tail']) for p,m in [('state','defect'),('state','source_error'),('heldout','defect'),('heldout','source_error')]] for i in range(4)]);im=a[0].imshow(-np.log10(v),vmin=0,vmax=4,cmap='Blues');a[0].set_xticks(range(4),['State defect','State fidelity','GR defect','GR fidelity'],rotation=20);a[0].set_yticks(range(4),['R0','R1','R2','R3']);a[0].set_title('A  Tested-reference minimum tails');f.colorbar(im,ax=a[0],label='−log10 conditional tail');s.extend({'panel':'A',**r} for r in rr)
    claims=read(ROOT/'completion/claim_matrix.csv');text(a[1],'B  Retained conditional authority','Strict GS-only geometry\nConstruction-held-out GR validation\nS1 aggregate + weakest support\nR0–R3 tested-reference nonreproduction\nPCA matched effect directions\n\nBiological replication: three embryos.');s.extend({'panel':'B',**r} for r in claims)
    closed=freeze()['closed_claims'];text(a[2],'C  Closed / outside authority','No reopening of:\ncount significance; composition enrichment;\nbridge enrichment; R4/R5; ownership; S2.\n\nHistorical lane bank: sensitivity, not primary.\nBaseline and biology: descriptive only.');s.extend({'panel':'C','claim':c,'status':'CLOSED_NOT_REOPENED'} for c in closed)
    text(a[3],'D  Exact/global boundary','No exact/global correspondence is established.\nNo strict exact-descent theorem is inferred\nfrom historical-object counterexamples.\nNo population generalization is asserted.\n\nAll retained draws are IID, not chains.\nNo additional rescue tranche is required.');s.append({'panel':'D','status':'OUTSIDE_AUTHORITY','claim':'strict exact descent/global equivalence/external confirmation'});finish(6,f,s)
    print('Six actual vector figures, PNG previews and source CSVs generated',flush=True)
if __name__=='__main__':main()

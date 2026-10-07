"""Final independent-calibration source, exact figure names and scientific Results order."""
import re,shutil
from io08 import *
from names08 import *
P7=ROOT.parent/'rmmo_paper3_acceptance07';P5=ROOT.parent/'rmmo_paper3_referee05';P6=ROOT.parent/'rmmo_paper3_revision06'
def comparisons(bank,panel,metric):
    return [r for r in read(ROOT/'statistical_authority/final_maxT_comparison.csv') if r['bank']==bank and r['panel']==panel and r['metric']==metric]
def figures():
    for old,name in FIGURES.items():
        if old in ('Figure2','FigureS3','FigureS4'):continue
        for ext in ('pdf','png'):
            path=ROOT/'publication/figures'/f'{name}.{ext}';path.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(P7/'figures'/f'{old}.{ext}',path)
        source=read(P7/'figures'/f'{old}_authority.csv')
        if old in ('Figure3','Figure6','FigureS2'):
            final=read(ROOT/'statistical_authority/final_component_statistics.csv')
            for row in source:
                if not row.get('maxT_tail'):continue
                bank='spatial_overlap/S1' if old=='Figure3' else 'retention/'+row['level']
                candidates=[x for x in final if x['bank']==bank and x['component']==row['component'] and (old=='Figure3' or x['panel']==row['panel'] and x['metric']==row['metric'])];assert len(candidates)==1
                for key in ('Z','maxT_tail','T_sum','T_min','sum_tail','min_tail'):row[key]=candidates[0][key]
                row['standardization']='independent_calibration'
        if old=='Figure5':
            final=read(ROOT/'statistical_authority/final_component_statistics.csv')
            for row in source:
                if not row.get('maxT_tail'):continue
                bank={'hash_primary':'acceptance07/representation/hash_s0_q25','pca32':'representation/pca32','pca64':'representation/pca64'}[row['variant']]
                stat=next(x for x in final if x['bank']==bank and x['panel']=='heldout' and x['metric']=='defect' and x['component']==row['component'])
                for key in ('Z','maxT_tail','T_sum','T_min','sum_tail','min_tail'):row[key]=stat[key]
                row['standardization']='independent_calibration'
        write(f'publication/figures/{name}_authority.csv',source)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'figure.dpi':120});f,axes=plt.subplots(2,2,figsize=(11.2,8));a=axes.ravel();source=[];f.suptitle('Figure 2. Construction-held-out molecular validation',fontsize=15)
    specs=[('A  Landmark-edge GR agreement','edge_agreement','heldout','edge_distance'),('B  Matched GR coherence','matched_coherence/q25','heldout','defect'),('C  Degree-controlled GR fidelity','source_fidelity/q25','heldout','source_error')]
    for ax,(title,bank,panel,metric) in zip(a,specs):
        rr=[r for r in read(P5/bank/'component_statistics.csv') if not r.get('panel') or r['panel']==panel and r['metric']==metric];cmp=comparisons(bank,panel,metric);assert len(rr)==len(cmp)==3
        combined=[]
        for r,c in zip(rr,cmp):
            assert r['component']==c['component'];stat=next(x for x in read(ROOT/'statistical_authority/final_component_statistics.csv') if x['bank']==bank and x['panel']==panel and x['metric']==metric and x['component']==r['component']);combined.append({**r,**{k:stat[k] for k in ('Z','maxT_tail','T_sum','T_min','sum_tail','min_tail')},'maxT_standardization':'independent_calibration','primitive_raw_tail_unchanged':True})
        x=np.arange(3);ax.bar(x-.16,[float(r['observed']) for r in rr],.32,label='Observed');ax.bar(x+.16,[float(r['null_mean']) for r in rr],.32,yerr=[float(r['null_sd']) for r in rr],capsize=3,label='Inferential null mean +/- draw SD');ax.set_xticks(x,[r['component'] for r in rr]);ax.set(title=title,ylabel='Molecular chord error/defect');ax.legend(fontsize=7);source.extend({'panel':title[0],**r,'plot_panel':title[0]} for r in combined)
        write(f'publication/figures/{FIGURES["Figure2"]}_{bank.split("/")[0]}.csv',combined)
    fidelity=comparisons('source_fidelity/q25','heldout','source_error');fam=next(r for r in read(ROOT/'statistical_authority/final_family_statistics.csv') if r['bank']=='source_fidelity/q25' and r['panel']=='heldout' and r['metric']=='source_error');a[3].axis('off');a[3].set_title('D  Unified calibration for the primary family',loc='left');a[3].text(0,.93,f"All three fidelity reductions positive.\nCalibrated sum tail = {float(fam['sum_tail']):.6g}\nCalibrated minimum tail = {float(fam['min_tail']):.6g}\nCalibrated O1 maxT = {float(fidelity[0]['calibrated_maxT_tail']):.6g}\nO1 remains individually unsupported.\n\nFamily support does not imply marginal support.\nDraw SD is not biological-replicate uncertainty.",va='top',fontsize=10,transform=a[3].transAxes);source.append({'panel':'D',**fam,'O1_maxT_tail':fidelity[0]['calibrated_maxT_tail']});f.tight_layout(rect=(0,0,1,.94));path=ROOT/'publication/figures';name=FIGURES['Figure2'];f.savefig(path/f'{name}.pdf');f.savefig(path/f'{name}.png',dpi=160);plt.close(f);write(f'publication/figures/{name}_authority.csv',source)
    scales=read(ROOT/'statistical_authority/final_scale_sensitivity.csv');f,axes=plt.subplots(2,2,figsize=(11.2,8));x=np.arange(3)
    for col,(bank,label) in enumerate([('matched_coherence','Matched GR coherence'),('source_fidelity','Degree-controlled GR fidelity')]):
        for component,_,_ in OVERLAPS:
            rr=[next(r for r in scales if r['bank']==bank+'/'+q and r['component']==component) for q in SCALES]
            axes[0,col].plot(x,[100*float(r['relative_effect']) for r in rr],marker='o',label=component)
        for stat in ('sum','min'):
            rr=[next(r for r in scales if r['bank']==bank+'/'+q and r['component']=='O1') for q in SCALES]
            axes[1,col].plot(x,[float(r[stat+'_tail']) for r in rr],marker='o',label=stat)
        axes[0,col].set(title=label,ylabel='Reduction (%)');axes[1,col].set(ylabel='Independently calibrated family tail',yscale='log')
        for ax in axes[:,col]:ax.set_xticks(x,list(SCALES));ax.legend(fontsize=8)
    f.suptitle('Figure S3. Frozen domain-scale sensitivity');f.tight_layout();name=FIGURES['FigureS3'];f.savefig(path/f'{name}.pdf');f.savefig(path/f'{name}.png',dpi=160);plt.close(f);write(f'publication/figures/{name}_authority.csv',scales)
    retention=read(ROOT/'statistical_authority/final_retention_effect_surface.csv');labels=[(panel,metric,component) for panel in ('state','heldout') for metric in ('defect','source_error') for component,_,_ in OVERLAPS]
    matrix=np.array([[float(next(r['Z'] for r in retention if (r['panel'],r['metric'],r['component'])==label and r['level']==level)) for level in ('R0','R1','R2','R3')] for label in labels]);f,ax=plt.subplots(figsize=(8.5,8));im=ax.imshow(matrix,aspect='auto',cmap='viridis');ax.set_xticks(range(4),['R0','R1','R2','R3']);ax.set_yticks(range(12),[(('GS' if p=='state' else 'GR')+' '+('D' if m=='defect' else 'E')+' '+c) for p,m,c in labels]);ax.set_title('Figure S4. Independently calibrated retention effects')
    for i in range(12):
        for j in range(4):ax.text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',color='white' if matrix[i,j]<matrix.max()*.6 else 'black',fontsize=8)
    f.colorbar(im,ax=ax,label='Independent-calibration observed Z');f.tight_layout();name=FIGURES['FigureS4'];f.savefig(path/f'{name}.pdf');f.savefig(path/f'{name}.png',dpi=160);plt.close(f);write(f'publication/figures/{name}_authority.csv',retention)
    write('publication/figures/name_registry.csv',[dict(previous=old+'.pdf',publication=name+'.pdf',source_data=name+'_authority.csv',updated_primary_calibration=old=='Figure2') for old,name in FIGURES.items()])
def sanitize(text):
    checks={
      'html_space_entity':r'&#x20;',
      'markdown_inside_textit':r'\\textit\{\*[^}]*\*\}',
      'doubled_math_parentheses':r'(?<!\\)\\\\[()]',
      'doubled_named_control':r'(?<!\\)\\\\(?:section|subsection|textit|textbf|begin|end|includegraphics|caption|Cref|cref|label|usepackage|documentclass)\b',
      'markdown_heading':r'(?m)^\s*#{1,6}\s',
      'markdown_code_fence':r'(?m)^\s*```',
      'markdown_link':r'\[[^\]\n]+\]\(https?://',
      'markdown_bold':r'\*\*[^*\n]+\*\*'}
    checks['markdown_inline_code']=r'(?<!`)`(?!`)[^`\n]+(?<!`)`(?!`)'
    matches={k:re.findall(v,text) for k,v in checks.items()};assert not any(matches.values()),matches
    return dict(actual_source_scanned=True,artifacts=matches,passed=True)
def manuscript():
    text=(P7/'manuscript_authority/Paper3_v9.tex').read_text();decision=json.loads((ROOT/'statistical_authority/final_calibration_decision.json').read_text());assert decision['all_reported_families_independently_calibrated'] and not decision['component_decisions_changed'] and not decision['family_decisions_changed']
    a=text.index('\\subsection{Strict construction and independent statistical closure}');b=text.index('\\subsection{The nine-lane construction',a);old=text[a:b];first=old[old.index('The final publication-primary'):old.index('The independent family-statistic')];table=old[old.index('\\begin{table}'):old.index('\\end{table}')+len('\\end{table}')]
    text=text[:a]+'\\subsection{Strict state-panel construction}\n\n'+first+'\\FloatBarrier\n'+text[b:]
    robust=r'''\FloatBarrier
\subsection{Statistical robustness under one independent-calibration convention}

Independent 1,024-draw same-law calibration banks supply fixed component means and sample SDs for every reported maxT family. Existing Revision06 calibration covers the 22 publication-primary families; seven additional calibration-only banks cover the 14 reported state/GR scale and representation sensitivity families. Across all 36 reported families, all 108 component maxT decisions and all 72 aggregate/minimum decisions are unchanged at the declared 0.05 threshold. No primitive statistic, null law, source geometry, effect size or inferential bank changed. Historical maxT values remain comparison authority only in the companion. The calibration plan, disjoint seed namespace and then the estimated constants were committed before the new outcomes and final ranks, respectively; no rescue was performed.

'''+table+r'''

The reported scope includes strict q10/q25/q50 matched coherence and degree-controlled source fidelity, hash-S0/PCA32/PCA64 q25 matched-defect controls, edge agreement, S1 and R0--R3. Target-PCA and companion-only cross-metric alternatives are not reported inferential families in this release and received no new calibration. The comparison companion retains their historical provenance separately.

'''
    marker='\\FloatBarrier\n\\section{Discussion}';assert marker in text;text=text.replace(marker,robust+marker,1)
    for old,new in FIGURES.items():text=text.replace('{'+old+'.pdf}','{'+new+'.pdf}')
    # The primary O1 number actually changes slightly even though its decision does not.
    row=comparisons('source_fidelity/q25','heldout','source_error')[0];text=text.replace('0.1220',f'{float(row["calibrated_maxT_tail"]):.4f}').replace('unsupported original O1 maxT component','unsupported independently calibrated O1 maxT component')
    a=text.index('\\subsection{Original component and maxT inferential standardization}');b=text.index('\\subsection{Joint family-draw construction}',a)
    methods=r'''\subsection{Independent calibration of component and family scores}

For every reported family, an independent same-law calibration bank of \(B_{\mathrm{cal}}=1024\) draws supplies fixed \(\mu_{i,\mathrm{cal}}\) and sample \(\sigma_{i,\mathrm{cal}}\), with divisor \(B_{\mathrm{cal}}-1\). Using the unchanged \(B=1024\) synchronized inferential draws, set \(s_i=-1\) for defects/errors and \(s_i=+1\) for overlap size:
\begin{align}
Z_{i,\mathrm{cal}}^{(b)}&=\frac{s_i(T_i^{(b)}-\mu_{i,\mathrm{cal}})}{\sigma_{i,\mathrm{cal}}+10^{-8}}, &
Z_{i,\mathrm{cal}}^{\mathrm{obs}}&=\frac{s_i(T_i^{\mathrm{obs}}-\mu_{i,\mathrm{cal}})}{\sigma_{i,\mathrm{cal}}+10^{-8}},\\
M_{\mathrm{cal}}^{(b)}&=\max_i Z_{i,\mathrm{cal}}^{(b)}, &
p_{i,\mathrm{maxT,cal}}&=\frac{1+\#\{b:M_{\mathrm{cal}}^{(b)}\ge Z_{i,\mathrm{cal}}^{\mathrm{obs}}\}}{1025}.
\end{align}
The maximum ranges over the three declared components of a family. Ties count as exceedances. The same calibrated component scores define \(T_\Sigma=\sum_i Z_{i,\mathrm{cal}}\) and \(T_{\min}=\min_i Z_{i,\mathrm{cal}}\); their synchronized inclusive plus-one tail ranks also use the original inferential draws. Primitive raw tails use the unchanged signed primitive statistics, without altered standardization.

Seven calibration-only banks were authorized to fill the missing reported q10/q50 and hash-S0/PCA32/PCA64 q25 families. The frozen namespace is \texttt{RMMO\_P3\_FINAL08\_CALIBRATION\_V1}; its 33,792 distinct seeds have no collision with the 122,880 prior audited seeds. Five matched-target banks each synchronize six orientation draws. The two degree-scale banks share the same 1,024 independent whole-graph realizations for each of the three pairs, as in the original cross-scale coupling; uniform stub permutations are rejected only for nonsimplicity. Matched-target sampling preserves each sourcewise target coarse-state/anchor-distance-decile histogram, uses distinct targets without replacement and uses direct uniform subset sampling for large requests. Every inferential bank remains byte-identical. No new inferential bank, null law, scale, representation, effect-size gate or rescue analysis was introduced. Independent calibration fixes normalization; the frozen inferential bank supplies the ranks.

\FloatBarrier
'''
    text=text[:a]+methods+text[b:]
    text=text.replace('All 22 publication families and all 44 aggregate/minimum support decisions remain unchanged.','All 36 reported families, all 108 component maxT decisions and all 72 aggregate/minimum support decisions remain unchanged.')
    text=text.replace('Primitive component tests and component maxT adjustment remain based on the original inferential law and are unchanged by the independent-calibration robustness check.','Primitive raw tails remain unchanged. Final component maxT, aggregate and minimum ranks all use independent-calibration standardization with the unchanged inferential draws.')
    text=text.replace('family tails and O1 maxT use original inferential-bank standardization','family tails and O1 maxT use independently fixed calibration constants')
    text=text.replace('uses existing banks, not new sampling or scale selection','uses unchanged inferential banks with independent calibration, not new inferential sampling or scale selection')
    text=text.replace('all 48 original-bank standardized retention effects','all 48 independently calibrated retention effects')
    text=text.replace('Original-bank aggregate and minimum tails','Independently calibrated aggregate and minimum tails')
    text=text.replace('All 48 original-bank standardized retention effects.','All 48 independently calibrated retention effects.')
    text=text.replace('Each cell is the existing signed observed Z using its original inferential null mean and sample draw SD plus','Each cell is the signed observed Z using its independent calibration mean and sample SD plus')
    text=text.replace('This is not calibrated-family Z or biological-replicate uncertainty.','This is calibrated-family Z, not biological-replicate uncertainty.')
    text=text.replace('Original-bank hash-S0 synchronized three-component family statistics. There is no new effect-size gate or independent calibration bank.','Independently calibrated hash-S0 synchronized three-component family statistics. There is no new effect-size gate; the inferential draws are unchanged.')
    final=read(ROOT/'statistical_authority/final_component_statistics.csv');families=read(ROOT/'statistical_authority/final_family_statistics.csv')
    for prefix,bank,metric in [('Matched GR D','matched_coherence','defect'),('Degree GR E','source_fidelity','source_error')]:
        for q in SCALES:
            rr=[r for r in final if r['bank']==bank+'/'+q and r['panel']=='heldout' and r['metric']==metric]
            replacement=prefix+' & '+q+' & '+' & '.join(f'{100*float(r["relative_effect"]):.2f}' for r in rr)+' & '+' & '.join(f'{float(rr[0][key]):.6f}' for key in ('sum_tail','min_tail','maxT_tail'))+r'\\'
            text=re.sub(r'(?m)^'+re.escape(prefix+' & '+q)+r' & .*$',lambda _:replacement,text)
    for panel,label in [('state','GS'),('heldout','GR')]:
        ff=next(r for r in families if r['bank']=='acceptance07/representation/hash_s0_q25' and r['panel']==panel)
        replacement=label+' & '+' & '.join(f'{float(ff[key]):.6f}' for key in ('T_sum','T_min','sum_tail','min_tail'))+r'\\'
        text=re.sub(r'(?m)^'+label+r' & [-\d.]+ & [-\d.]+ & [-\d.]+ & [-\d.]+\\\\$',lambda _:replacement,text)
    text=text.replace('A separate 1,024-draw calibration bank fixes component standardization','Existing separate 1,024-draw calibration banks fix component standardization')
    text=text.replace(r'\texttt{RMMO\_P3\_FINAL08\_CALIBRATION\_V1}',r'\texttt{RMMO\_}\allowbreak\texttt{P3\_}\allowbreak\texttt{FINAL08\_}\allowbreak\texttt{CALIBRATION\_}\allowbreak\texttt{V1}')
    text=text.replace('coarse-state/anchor-distance-decile','coarse state and anchor-distance decile')
    text=text.replace('hash-S0/PCA32/PCA64','hash-S0, PCA32 and PCA64')
    text=text.replace("`searchsorted(side='right')`",r"\texttt{searchsorted(side='right')}").replace('`scipy.sparse.linalg.svds`',r'\texttt{scipy.sparse.linalg.svds}')
    # Statements about original-bank sensitivity statistics remain true and explicit.
    receipt=sanitize(text);md(f'publication/manuscript/{MANUSCRIPT}.tex',text);js('publication/source_sanitization.json',{**receipt,'tex_sha256':sha(ROOT/f'publication/manuscript/{MANUSCRIPT}.tex'),'all_reported_families_calibrated':True});print('Final exact-name manuscript prepared with one calibrated convention',flush=True)
def main():figures();manuscript()
if __name__=='__main__':main()

"""Read-only final manuscript/table/figure validation, with optional private parent check."""
import argparse,re,subprocess
from pypdf import PdfReader
from PIL import Image
from io08 import *
from names08 import *
from presentation08 import sanitize
def validate(private=False):
    directory=ROOT/'publication';source=directory/'manuscript'/f'{MANUSCRIPT}.tex';pdf=source.with_suffix('.pdf');text=source.read_text();sanitize(text)
    assert not re.search(r'(?<!\\)\\\\[\[\]](?![\d.]+(?:em|pt|ex|cm|mm)\])',text)
    assert not any(x in text for x in ('original-bank standardiz','original inferential-bank standardization','Provisional Final08','coverage unresolved','maxT uses original','no separate calibration bank'))
    sections=['Introduction','Results','Discussion','Scope and limitations','Methods','Data availability'];assert [text.index('\\section{'+s+'}') for s in sections]==sorted(text.index('\\section{'+s+'}') for s in sections)
    for needle in ('\\documentclass[11pt]','margin=1in','\\setstretch{1.08}','\\hyphenpenalty=10000','\\captionsetup{font=small,labelfont=bf,textfont=it}','\\usepackage[nameinlink,capitalise,noabbrev]{cleveref}','\\vspace*{0.17\\textheight}'):
        assert needle in text,needle
    assert text.count('\\begin{figure}[H]')==10 and text.count('\\FloatBarrier')>20
    refs=re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',text);assert set(refs)=={n+'.pdf' for n in FIGURES.values()} and len(refs)==10
    log=source.with_suffix('.log').read_text();assert not any(x in log for x in ('Overfull \\hbox','Overfull \\vbox','undefined references','Float too large','Citation `'))
    pages=PdfReader(pdf).pages;assert len(pages)>20 and '??' not in ' '.join(x.extract_text() or '' for x in pages)
    final=read(ROOT/'statistical_authority/final_component_statistics.csv')
    for old,name in FIGURES.items():
        assert len(PdfReader(directory/'figures'/f'{name}.pdf').pages)==1 and Image.open(directory/'figures'/f'{name}.png').width>1000
        rr=read(directory/'figures'/f'{name}_authority.csv');assert rr
        for row in rr:
            if not row.get('maxT_tail'):continue
            assert (row.get('standardization') or row.get('maxT_standardization'))=='independent_calibration',old
            # Every figure-facing component statistic uses the final numeric authority.
            panel=row.get('panel');metric=row.get('metric');bank=row.get('bank')
            if old=='Figure2':bank={'A':'edge_agreement','B':'matched_coherence/q25','C':'source_fidelity/q25'}[row['plot_panel']];panel='heldout';metric={'A':'edge_distance','B':'defect','C':'source_error'}[row['plot_panel']]
            elif old=='Figure3':bank='spatial_overlap/S1';panel='spatial';metric='overlap_size'
            elif old in ('Figure6','FigureS2'):bank='retention/'+row['level']
            elif old=='Figure5':bank={'hash_primary':'acceptance07/representation/hash_s0_q25','pca32':'representation/pca32','pca64':'representation/pca64'}[row['variant']];panel='heldout';metric='defect'
            stat=next(x for x in final if x['bank']==bank and x['panel']==panel and x['metric']==metric and x['component']==row['component'])
            for key in ('Z','maxT_tail','T_sum','T_min','sum_tail','min_tail'):assert float(row[key])==float(stat[key]),(old,key)
    assert len(read(directory/'figures'/f'{FIGURES["FigureS3"]}_authority.csv'))==18
    assert len(read(directory/'figures'/f'{FIGURES["FigureS4"]}_authority.csv'))==48
    parent_count=None
    if private:
        from verify08 import parents
        parent_count=parents()
        for receipt,files in [('missing_calibration_commit.json',['freeze/missing_calibration_plan.json','banks08.py']),('constants_commit_receipt.json',['statistical_authority/component_calibration_final.csv'])]:
            commit=json.loads((ROOT/'freeze'/receipt).read_text())['reachable_commit']
            for rel in files:
                content=subprocess.check_output(['git','show',commit+':research/rmmo_paper3_final08/'+rel],cwd=ROOT.parent.parent);assert hashlib.sha256(content).hexdigest()==sha(ROOT/rel)
    js('completion/publication_validation.json',dict(passed=True,manuscript_pages=len(pages),figures=10,all_figure_sources_use_final_calibration=True,source_artifacts=0,overfull_boxes=0,undefined_references=0,methods_after_scope=True,exact_series_shell_preserved=True,parent_artifacts_unchanged=parent_count,pdf_sha256=sha(pdf),tex_sha256=sha(source)))
    print('PASS publication tables, exact names, sanitized source and PDF:',len(pages),'pages',flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--private',action='store_true');validate(parser.parse_args().private)

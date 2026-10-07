"""Original Final08 figure/statistical comparisons, without an older manuscript PDF."""
import csv
import re
import sys
from pathlib import Path
from pypdf import PdfReader
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'reproducibility/rmmo_paper3_final08'
sys.path.insert(0, str(BASE))
from names08 import FIGURES

def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))

def main():
    final = rows(BASE / 'statistical_authority/final_component_statistics.csv')
    refs = []
    for rel, expected in [('manuscript/RMMO_Paper3_Manuscript_preDOI.tex', 6),
                          ('supplement/SupplementaryEvidence.tex', 4)]:
        text = (ROOT / rel).read_text()
        found = re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', text)
        assert len(found) == expected
        refs.extend(found)
    assert len(refs) == 10 and set(refs) == {n+'.pdf' for n in FIGURES.values()}
    checked = 0
    for old, name in FIGURES.items():
        assert len(PdfReader(ROOT / 'figures' / (name+'.pdf')).pages) == 1
        with Image.open(ROOT / 'figures' / (name+'.png')) as image:
            image.verify()
        with Image.open(ROOT / 'figures' / (name+'.png')) as image:
            assert image.width > 1000
        rr = rows(ROOT / 'figures' / (name+'_authority.csv'))
        assert rr
        for row in rr:
            if not row.get('maxT_tail'):
                continue
            assert (row.get('standardization') or row.get('maxT_standardization')) == 'independent_calibration', old
            panel, metric, bank = row.get('panel'), row.get('metric'), row.get('bank')
            if old == 'Figure2':
                bank = {'A':'edge_agreement', 'B':'matched_coherence/q25', 'C':'source_fidelity/q25'}[row['plot_panel']]
                panel = 'heldout'
                metric = {'A':'edge_distance', 'B':'defect', 'C':'source_error'}[row['plot_panel']]
            elif old == 'Figure3':
                bank, panel, metric = 'spatial_overlap/S1', 'spatial', 'overlap_size'
            elif old in ('Figure6', 'FigureS2'):
                bank = 'retention/' + row['level']
            elif old == 'Figure5':
                bank = {'hash_primary':'acceptance07/representation/hash_s0_q25', 'pca32':'representation/pca32', 'pca64':'representation/pca64'}[row['variant']]
                panel, metric = 'heldout', 'defect'
            stat = next(x for x in final if x['bank'] == bank and x['panel'] == panel
                        and x['metric'] == metric and x['component'] == row['component'])
            for key in ('Z', 'maxT_tail', 'T_sum', 'T_min', 'sum_tail', 'min_tail'):
                assert float(row[key]) == float(stat[key]), (old, key)
                checked += 1
    assert len(rows(ROOT / 'figures' / (FIGURES['FigureS3']+'_authority.csv'))) == 18
    assert len(rows(ROOT / 'figures' / (FIGURES['FigureS4']+'_authority.csv'))) == 48
    print(f'PASS ten frozen vector PDFs/PNGs, six main/four supplement references and {checked} original Final08 figure-statistic comparisons')

if __name__ == '__main__':
    main()

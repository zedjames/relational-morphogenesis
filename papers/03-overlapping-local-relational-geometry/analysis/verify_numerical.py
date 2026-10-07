"""Full Final08 numerical gate, independent of scholarly PDF/layout availability."""
import argparse
import json
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path
from verify_artifacts import ROOT, verify

def run(script, args, env):
    subprocess.run([sys.executable, str(script)] + args, check=True, env=env)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--receipt', type=Path, help='Optional successful-run receipt outside frozen authorities')
    args = parser.parse_args()
    export = verify()
    base = ROOT / 'reproducibility'
    with tempfile.TemporaryDirectory(prefix='rmmo-p3-replay-') as runtime:
        env = dict(os.environ, RMMO_P3_INPUT_ROOT=str(base / 'data_authority'),
                   RMMO_P3_CACHE=str(Path(runtime) / 'cache'),
                   PYTHONPYCACHEPREFIX=str(Path(runtime) / 'pycache'),
                   OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', OMP_NUM_THREADS='1')
        if platform.system() == 'Linux' and platform.machine() == 'x86_64':
            # Pin execution kernels, not scientific algorithms or tolerances.
            flags = Path('/proc/cpuinfo').read_text().split()
            if not {'avx2', 'fma'}.issubset(flags):
                raise RuntimeError('Frozen Linux replay requires AVX2/FMA for the pinned Haswell OpenBLAS kernel')
            env['OPENBLAS_CORETYPE'] = 'HASWELL'
            print('Runtime: Linux x86_64 / single-thread OpenBLAS HASWELL', flush=True)
        run(ROOT / 'analysis/verify_binary_release.py', [], env)
        for action in ('scalar_ranks', 'matched_replay', 'degree_replay'):
            run(base / 'rmmo_paper3_final08/verify_calibration08.py', [action], env)
        run(ROOT / 'analysis/verify_geometry.py', [], env)
        run(ROOT / 'analysis/verify_figures.py', [], env)
    ledger = json.loads((base / 'rmmo_paper3_final08/statistical_authority/final_calibration_decision.json').read_text())
    assert ledger['status'] == 'CALIBRATION_PRESERVED' and ledger['publication_freeze_complete']
    assert ledger['component_decisions_changed'] == ledger['family_decisions_changed'] == []
    assert ledger['old_inferential_banks_changed'] == 0 and ledger['no_rescue']
    verify()
    receipt = {'passed': True, 'scientific_files': export['exported_files'],
               'scientific_bytes': export['exported_bytes'], 'independent_calibration_banks': 7,
               'exact_graph_draw_computations': 33792, 'components': 108,
               'families': 36, 'aggregate_minimum_events': 72,
               'component_decisions_changed': 0, 'family_decisions_changed': 0,
               'inferential_banks_changed': 0, 'figures': 10,
               'exact_composition_routes': 3, 'exact_bridge_overlaps': 3,
               'private_checkout_required': False, 'raw_acquisition_required': False,
               'document_layout_validation': 'separate'}
    if args.receipt:
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print('PASS full standalone Final08 numerical reproducibility: ' + json.dumps(receipt, sort_keys=True))

if __name__ == '__main__':
    main()

"""Read-only first-draw diagnostics; never a substitute for the full replay gate."""
import hashlib
import json
import platform
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'reproducibility/rmmo_paper3_final08'))
from banks08 import Sampler, context, rng, p, read

def digest(value):
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()

def main():
    bank = 'acceptance07/representation/hash_s0_q25'
    rows, z, ix, cores, geometry, overlaps, matrices = p.arrays('hash_primary')
    records = {r['orientation']: r for r in read(
        ROOT / 'reproducibility/rmmo_paper3_final08/independent_calibration_bank' /
        bank / 'graph_audit.csv') if int(r['draw']) == 0}
    output = {'python': platform.python_version(), 'os': platform.system(),
              'os_release': platform.release(), 'machine': platform.machine(),
              'numpy': np.__version__,
              'blas': np.show_config(mode='dicts')['Build Dependencies']['blas']['name'],
              'state_sha256': digest(z), 'bank': bank, 'orientations': []}
    for name, edges in geometry.items():
        sampler = Sampler(2, rows, z, ix, edges, 'q25')
        strata, lookup, observed = context(sampler)
        sampler.draw(rng(f'matched|{bank}|0|{name}'), {})
        output['orientations'].append({
            'orientation': name, 'strata': strata,
            'stratum_lookup_sha256': digest(lookup),
            'observed_membership_sha256': digest(sampler.observed),
            'deciles_sha256': digest(sampler.dec),
            'expected_graph_hash': records[name]['graph_hash'],
            'actual_graph_hash': sampler.last_graph_hash,
            'exact_graph_match': sampler.last_graph_hash == records[name]['graph_hash']})
    print(json.dumps(output, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()

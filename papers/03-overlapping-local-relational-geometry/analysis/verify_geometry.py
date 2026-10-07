"""Read-only replay of frozen strict geometry, composition and existential bridges."""
import csv
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'reproducibility'
sys.path.insert(0, str(BASE / 'rmmo_paper3_referee05'))
import common05 as p

def main():
    r, z, ix, cores, geometry, overlaps, matrices = p.arrays()
    assert len(r) == 3704 and [len(x) for x in ix] == [1032, 1298, 1374]
    assert [len(cores[k][0]) for k in p.PAIRS] == [59, 17, 44]
    lanes = p.read(p.ROOT / 'strict_lanes/consensus_edges.csv.gz')
    assert len(lanes) == 120 and all(7 <= int(row['lane_support']) <= 9 for row in lanes)
    saved = np.load(p.ROOT / 'geometry/stageA_arrays.npz')
    relations, fibers = {}, {}
    for name, e in geometry.items():
        for q in p.SCALES:
            np.testing.assert_array_equal(e[q]['nodes'], saved[name+'_'+q+'_nodes'])
            np.testing.assert_array_equal(e[q]['active'], saved[name+'_'+q+'_active'])
            np.testing.assert_allclose(e[q]['radius'], float(saved[name+'_'+q+'_radius']), rtol=2e-11, atol=2e-12)
            targets, membership = p.fiber_matrix(e, q)
            np.testing.assert_array_equal(targets, saved[name+'_'+q+'_fiber_targets'])
            np.testing.assert_array_equal(membership, saved[name+'_'+q+'_fiber_membership'])
            if q == 'q25':
                fibers[name] = {int(n): set(map(int, targets[membership[j]])) for j, n in enumerate(e[q]['nodes'])}
                relations[name] = {(n, t) for n, ts in fibers[name].items() for t in ts}
    for (name, q), cells in overlaps.items():
        np.testing.assert_array_equal(cells, saved[name+'_'+q])
    authority = BASE / 'rmmo_paper3_revision06/methods_authority'
    for row in p.read(authority / 'composition_exact_counts.csv'):
        a, b, c = row['route'].split('->')
        overlap = {'1': 'O1', '2': 'O2', '3': 'O3'}[a]
        sources = set(map(int, overlaps[overlap, 'q25']))
        direct = {(x, t) for x, t in relations['L'+a+c] if x in sources}
        composed = {(x, t) for x in sources for y in fibers['L'+a+b][x]
                    for t in fibers['L'+b+c].get(y, ())}
        actual = {'direct': len(direct), 'composed': len(composed),
                  'intersection': len(direct & composed), 'direct_only': len(direct-composed),
                  'composed_only': len(composed-direct), 'exact_equality': str(direct == composed)}
        assert all(str(value) == row[key] for key, value in actual.items()), (row['route'], actual)
    cellwise = {(row['overlap'], int(row['source'])): row for row in p.read(authority / 'refinement_cellwise.csv')}
    for row in p.read(authority / 'refinement_exact_counts.csv'):
        name = row['overlap']
        _, first, second = next(s for s in p.OVERLAPS if s[0] == name)
        b, c = first[-1], second[-1]
        pair = 'P' + ''.join(sorted((b, c)))
        left, right = cores[pair]
        bridges = set(zip(map(int, left), map(int, right))) if b < c else set(zip(map(int, right), map(int, left)))
        satisfying = coverage = 0
        for x0 in overlaps[name, 'q25']:
            x = int(x0)
            fa, fb = fibers[first][x], fibers[second][x]
            found = bridges & {(a, b0) for a in fa for b0 in fb}
            full = {a for a, _ in found} == fa and {b0 for _, b0 in found} == fb
            recorded = cellwise[name, x]
            assert int(recorded['bridge_pairs']) == len(found)
            assert recorded['historical_exact_criterion'] == str(bool(found))
            assert recorded['full_partner_coverage'] == str(full)
            satisfying += bool(found)
            coverage += full
        tested = len(overlaps[name, 'q25'])
        assert tested == int(row['tested']) and satisfying == int(row['satisfying'])
        assert tested-satisfying == int(row['failing']) and coverage == int(row['full_partner_coverage'])
        assert satisfying/tested == float(row['fraction']) and row['complete'] == str(satisfying == tested)
    print('PASS frozen 3704-cell / nine-lane / 59-17-44 geometry; three exact compositions and all 746 cellwise bridge diagnostics')

if __name__ == '__main__':
    main()

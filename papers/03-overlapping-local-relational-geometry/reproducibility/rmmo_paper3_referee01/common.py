"""Owned deterministic numerics; historical inputs are strictly read-only."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.stats import binom

ROOT = Path(__file__).resolve().parent
RESEARCH = Path(os.environ.get('RMMO_P3_INPUT_ROOT', str(ROOT.parent/'data_authority')))
CACHE = Path(os.environ.get('RMMO_P3_CACHE', str(Path(__import__('tempfile').gettempdir())/'rmmo_p3_referee01_cache')))
B = 1024
REPS = ('E7.5-R1', 'E7.5-R2', 'E7.5-R3')
PAIRS = ('P12', 'P13', 'P23')
ENDS = {'P12': (0, 1), 'P13': (0, 2), 'P23': (1, 2)}
SCALES = {'q10': .1, 'q25': .25, 'q50': .5}
OVERLAPS = (('O1', 'L12', 'L13'), ('O2', 'L21', 'L23'), ('O3', 'L31', 'L32'))
DELTA = .29931360483169556


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''): h.update(chunk)
    return h.hexdigest()


def seed(label):
    return int.from_bytes(hashlib.sha256(('RMMO_P3_REFEREE01_V1|' + label).encode()).digest()[:8], 'big')


def rng(label): return np.random.default_rng(seed(label))


def read(path):
    path = Path(path)
    with (gzip.open if path.suffix == '.gz' else open)(path, 'rt', newline='') as f:
        return list(csv.DictReader(f))


def write(path, rows, fields=None):
    path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    assert fields, path
    with (gzip.open if path.suffix == '.gz' else open)(path, 'wt', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def js(path, value):
    path = ROOT / path; path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def check_freeze():
    receipt = json.loads((ROOT / 'freeze/freeze_receipt.json').read_text())
    for path, expected in receipt['file_sha256'].items():
        assert sha(ROOT / path) == expected, ('plan changed', path)
    return receipt


def parent(t, path): return RESEARCH / f'relational_morphogenesis{t}' / path


def unit(a):
    a = np.asarray(a, dtype=np.float64)
    n = np.linalg.norm(a, axis=-1, keepdims=True)
    return a / np.where(n == 0, 1, n)


def chord(a, b):
    with np.errstate(all='ignore'):
        d = np.sqrt(np.maximum(0., 2. - 2. * (a @ b.T)))
    assert np.isfinite(d).all()
    return d


def load_state():
    rows = read(parent(12, 'authority/transcriptome_lineage_join.csv.gz'))
    with np.load(parent(12, 'embedding/state_projection.npz')) as f: zs = unit(f['ZS'])
    assert len(rows) == 3704
    return rows, zs


def load_target():
    with np.load(parent(12, 'embedding/realization_projection.npz')) as f: return f['ZR'].astype(np.float64)


def indices(rows): return [np.flatnonzero([r['replicate'] == rep for r in rows]) for rep in REPS]


def npz_parent(t, path):
    with np.load(parent(t, path)) as f: return {k: f[k].copy() for k in f.files}


def cores():
    bank = npz_parent(15, 'pair_native/pair_native_cores.npz')
    return {p: (bank[p + '_left'], bank[p + '_right']) for p in PAIRS}


def pair_distances(z, ix):
    ds, breaks, bins = {}, {}, {}
    for p, (a, b) in ENDS.items():
        ds[p] = chord(z[ix[a]], z[ix[b]])
        breaks[p] = np.quantile(ds[p], np.arange(.1, 1., .1))
        bins[p] = np.searchsorted(breaks[p], ds[p], side='right').astype(np.uint8)
    return ds, breaks, bins


def geometry(z, ix, core):
    """All source geometry constructed without target panel, annotations/history."""
    out = {}
    for p, (a, b) in ENDS.items():
        left, right = core[p]
        for source, target, edge_s, edge_t in [(a, b, left, right), (b, a, right, left)]:
            name = f'L{source+1}{target+1}'
            anchors = np.unique(edge_s)
            assert len(anchors), (p, name)
            d = chord(z[ix[source]], z[anchors])
            rho = d.min(1)
            entry = {'source': source, 'target': target, 'anchors': anchors, 'edge_s': edge_s, 'edge_t': edge_t, 'rho': rho, 'all_source': ix[source]}
            for scale, q in SCALES.items():
                radius = float(np.quantile(rho, q))
                mask = rho <= radius
                entry[scale] = {'nodes': ix[source][mask], 'radius': radius, 'active': d[mask] <= radius}
            out[name] = entry
    overlaps = {(o, s): np.intersect1d(out[a][s]['nodes'], out[b][s]['nodes']) for o, a, b in OVERLAPS for s in SCALES}
    return out, overlaps


def fiber_matrix(entry, scale, target_edges=None, source_edges=None):
    """Rows on fixed source domain, columns on exact target endpoint set."""
    te = entry['edge_t'] if target_edges is None else np.asarray(target_edges)
    targets = np.unique(te)
    adj = np.zeros((len(entry['anchors']), len(targets)), dtype=np.int32)
    se = entry['edge_s'] if source_edges is None else np.asarray(source_edges)
    adj[np.searchsorted(entry['anchors'], se), np.searchsorted(targets, te)] = 1
    # int32 prevents boolean/wrapping-count multiplication from changing fibers.
    membership = (entry[scale]['active'].astype(np.int32) @ adj) > 0
    assert membership.any(1).all(), 'nonempty is a construction invariant'
    return targets, membership


def means(matrix, targets, membership, nodes, domain):
    loc = np.searchsorted(domain, nodes)
    weights = membership[loc].astype(np.float64)
    with np.errstate(all='ignore'):
        values = weights @ matrix[targets]
    sizes = weights.sum(1)
    assert (sizes > 0).all() and np.isfinite(values).all()
    return values / sizes[:, None]


def score(first, second, source, normalized):
    if normalized: first, second = unit(first), unit(second)
    defect = np.linalg.norm(first - second, axis=1)
    error = (np.linalg.norm(source-first, axis=1) + np.linalg.norm(source-second, axis=1))/2
    assert np.isfinite(defect).all() and np.isfinite(error).all()
    return {'defect': float(defect.mean()), 'source_error': float(error.mean()), 'tolerance_fraction': float(np.mean(defect <= DELTA))}


def null_summary(observed, values, lower=False):
    values = np.asarray(values, dtype=float)
    assert len(values) == B and np.isfinite(values).all()
    ordered = np.sort(values)
    lo = max(0, int(binom.ppf(.025, B, .95))-1)
    hi = min(B-1, int(binom.ppf(.975, B, .95)))
    return {'observed': float(observed), 'null_mean': float(values.mean()), 'null_sd': float(values.std(ddof=1)), 'null95': float(np.quantile(values, .95)), 'null95_lo': float(ordered[lo]), 'null95_hi': float(ordered[hi]), 'conditional_tail': float((1 + np.sum(values <= observed if lower else values >= observed))/(B+1)), 'reduction': float((values.mean()-observed)/values.mean()) if values.mean() else 0.}


def family(observed, null, lower=False):
    """Rows synchronized replicates; columns family members. Finite reference tails."""
    null = np.asarray(null, dtype=float); obs = np.asarray(observed, dtype=float)
    assert null.shape[0] == B and null.shape[1] == len(obs)
    sign = -1 if lower else 1
    mu, sd = null.mean(0), null.std(0, ddof=1)+1e-8
    zz, zo = sign*(null-mu)/sd, sign*(obs-mu)/sd
    minimum, maximum = zz.min(1), zz.max(1)
    return {'all_three_minT': float(zo.min()), 'all_three_tail': float((1+np.sum(minimum >= zo.min()))/(B+1)), 'maxT_adjusted_tails': [float((1+np.sum(maximum >= x))/(B+1)) for x in zo], 'all_improved': bool(np.all(sign*(obs-mu) > 1e-8)), 'reference_authority': 'conditional finite-algorithm reference; not biological replication'}

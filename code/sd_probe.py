"""Numerical evidence for Conjecture 9: xi <= -2 psi for stochastically decreasing copulas.

Searches the checkerboard copulas on an n x n grid whose conditional distribution
functions are nondecreasing in the predictor (the stochastically decreasing
checkerboards) for large values of xi + lam * psi. For such copulas xi and psi are
evaluated exactly (cell densities are constant, so conditional distribution
functions are piecewise linear). xi + lam*psi is convex in the cell masses, so its
maximum over the polytope is attained at a vertex; the search hops between vertices
by repeated linear programming from seeded random starts. This is a local search:
it supplies evidence, not a proof.

Run from the package root:  python3 code/sd_probe.py --output results/sd_probe.json
Requires NumPy and SciPy (HiGHS). Deterministic for a fixed seed.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def setup(n):
    idx = np.arange(n)
    dist = abs(idx[:, None]-idx[None, :]).astype(float)
    dist[idx, idx] = 1/3
    a_eq, b_eq = [], []
    for i in range(n):
        a = np.zeros((n, n)); a[i, :] = 1; a_eq.append(a.ravel()); b_eq.append(1/n)
    for j in range(n-1):
        a = np.zeros((n, n)); a[:, j] = 1; a_eq.append(a.ravel()); b_eq.append(1/n)
    a_ub = []
    for i in range(n-1):
        for j in range(n-1):
            a = np.zeros((n, n)); a[i, :j+1] = 1; a[i+1, :j+1] = -1; a_ub.append(a.ravel())
    return dist, np.array(a_eq), np.array(b_eq), np.array(a_ub)


def stats(P, dist):
    n = len(P)
    H = n*np.cumsum(P, axis=1)
    L = np.c_[np.zeros(n), H[:, :-1]]
    xi = 6*np.sum((L*L+L*H+H*H)/3)/n/n - 2
    psi = 1 - 3*np.sum(dist*P)/n
    return xi, psi


def gradient(P, lam, dist):
    n = len(P)
    H = n*np.cumsum(P, axis=1)
    L = np.c_[np.zeros(n), H[:, :-1]]
    g = (2*H+L)/3
    g[:, :-1] += (2*L[:, 1:]+H[:, 1:])/3
    g *= 6/n/n
    return n*np.cumsum(g[:, ::-1], axis=1)[:, ::-1] - lam*3*dist/n


def search(n, lam, starts, rng):
    dist, a_eq, b_eq, a_ub = setup(n)
    b_ub = np.zeros(len(a_ub))

    def lp(c):
        r = linprog(c, A_ub=a_ub, b_ub=b_ub, A_eq=a_eq, b_eq=b_eq, bounds=(0, None), method='highs')
        if r.status != 0:
            raise ArithmeticError('linear program failed')
        return r.x.reshape(n, n)

    best = None
    for s in range(starts+1):
        # start 0 is independence (all cells equal); the rest are random vertices
        P = np.full((n, n), 1/(n*n)) if s == 0 else lp(rng.standard_normal(n*n))
        for _ in range(500):
            Q = lp(-gradient(P, lam, dist).ravel())
            fP = np.dot(stats(P, dist), [1, lam]); fQ = np.dot(stats(Q, dist), [1, lam])
            if fQ <= fP + 1e-13:
                break
            P = Q
        # feasibility of the returned vertex
        H = np.cumsum(P, axis=1)
        if np.min(np.diff(H, axis=0)) < -1e-7 or abs(P.sum(axis=0)-1/n).max() > 1e-7 or P.min() < -1e-9:
            raise ArithmeticError('returned point is not a stochastically decreasing checkerboard')
        xi, psi = stats(P, dist)
        if best is None or xi+lam*psi > best['value']:
            best = {'value': xi+lam*psi, 'xi': xi, 'psi': psi}
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grids', nargs='+', type=int, default=[10, 16, 24])
    parser.add_argument('--lams', nargs='+', type=float, default=[1.5, 1.9, 2.0, 2.1, 3.0])
    parser.add_argument('--starts', type=int, default=40)
    parser.add_argument('--seed', type=int, default=20261124)
    parser.add_argument('--output', default='results/sd_probe.json')
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    rows = []
    for n in args.grids:
        for lam in args.lams:
            b = search(n, lam, args.starts, rng)
            b.update({'n': n, 'lam': lam, 'xi_plus_2psi_at_best': b['xi']+2*b['psi']})
            rows.append(b)
            print(f"n={n} lam={lam}: max xi+lam*psi = {b['value']:.6f} at (xi,psi)=({b['xi']:.4f},{b['psi']:.4f})", flush=True)
    # psi <= 0 for SD copulas, so xi + lam*psi <= xi + 2*psi for lam >= 2: lam = 2 is the decisive case.
    worst = max(r['value'] for r in rows if r['lam'] >= 2)
    out = {'rows': rows, 'max_over_lam_ge_2': worst,
           'note': 'local vertex search; evidence for Conjecture 9, not a proof'}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2))
    print('max of xi + lam*psi found for lam >= 2:', f'{worst:.3e}')


if __name__ == '__main__':
    main()

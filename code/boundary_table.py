"""Tables and constants quoted in the manuscript.

Evaluates the boundary point (x_k, y_k) of Theorem 1 by the closed form of
Proposition 3 when k >= 1 and by the finite-mode solver (finite_modes.py)
when k <= 1, cross-checks the two at k = 1 and against the midpoint
Nystrom solver of verify_boundary.py, and writes
  * the least attainable xi at selected footrule values, with Rockel's
    Jensen bound and the Cauchy-Schwarz bound psi^2 for comparison;
  * the values of the stochastically decreasing family B_l (Proposition 8);
  * the constants quoted in the text.
Run from the package root:  python3 code/boundary_table.py --output results/table.json
Requires NumPy and SciPy. Every check raises on failure (safe under python -O).
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from finite_modes import boundary  # noqa: E402
from verify_boundary import curve  # noqa: E402


def closed_form(k):
    d = math.sqrt(2/k)
    D = math.cosh(d/2)
    x = 0.5 - 3/(8*D*D)*(math.cosh(d) - math.sinh(d)/d)
    y = -0.5 + 3*d*d/(32*D*D)*(math.sinh(d)/d - 1)
    return x, y


def point(k):
    if k >= 1:
        return closed_form(k)
    return boundary(k)


def k_for_psi(psi):
    if not -0.5 < psi < 0:
        raise ValueError(psi)
    return brentq(lambda k: point(k)[1] - psi, 0.02, 1e4, xtol=1e-15, rtol=1e-15, maxiter=500)


def rockel_jensen(psi):
    """Rockel's lower bound xi(C_mu) at footrule psi (JCAM 2026, Prop. 2 and Thm 3.2)."""
    def psi_of(v):
        return -2*v*v + 6*v - 5 + 1/v

    v = brentq(lambda v: psi_of(v) - psi, 0.5, 1.0, xtol=1e-15)
    return -4*v*v + 20*v - 17 + 2/v - 1/(v*v) - 12*math.log(v)


def block_family(l, n=4000):
    """xi and psi of B_l: independence on the central square of side l, V = 1-U elsewhere.

    Evaluated by direct quadrature of the conditional distribution function,
    independently of the closed forms 1-l^2 and -(1-l^2)/2.
    """
    a, b = (1-l)/2, (1+l)/2
    g = (np.arange(n)+0.5)/n
    U, V = np.meshgrid(g, g, indexing='ij')
    central = (U > a) & (U < b)
    h = np.where(central, np.clip((V-a)/l, 0, 1) if l > 0 else 0, (V >= 1-U).astype(float))
    xi = 6*np.mean((h-V)**2)
    # psi = 6 int C(t,t) dt - 2 with C(t,t) = int_0^t h(u,t) du, by the same quadrature.
    diag = np.mean(np.where(U <= V, h, 0.0), axis=0)
    psi = 6*np.mean(diag) - 2
    # SD check: h nondecreasing in u for every v.
    if np.min(np.diff(h, axis=0)) < -1e-12:
        raise ArithmeticError('B_l is not stochastically decreasing')
    return xi, psi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/table.json')
    args = parser.parse_args()
    out = {}

    # Cross-checks of the two evaluation methods.
    fm1, cf1 = boundary(1.0), closed_form(1.0)
    diff = max(abs(fm1[0]-cf1[0]), abs(fm1[1]-cf1[1]))
    if diff > 1e-10:
        raise ArithmeticError(f'finite modes and closed form disagree at k=1: {diff}')
    nystrom = {}
    for k in (0.125, 0.25, 0.5, 1, 2, 4, 8):
        exact = point(k)
        approx = curve(k, 4096)
        err = max(abs(exact[0]-approx[0]), abs(exact[1]-approx[1]))
        if err > 1e-6:
            raise ArithmeticError(f'midpoint solver disagrees at k={k}: {err}')
        nystrom[str(k)] = err
    out['cross_checks'] = {'finite_modes_vs_closed_form_k1': diff,
                           'midpoint_4096_vs_exact_max_abs': nystrom}

    # Constants quoted in the text.
    x1, y1 = point(1.0)
    xh, yh = point(0.5)
    out['constants'] = {
        'k1_xi': x1, 'k1_psi': y1, 'k1_xi_plus_2psi': x1 + 2*y1,
        'khalf_xi': xh, 'khalf_psi': yh, 'khalf_xi_plus_psi': xh + yh,
        'k0p3_xi': point(0.3)[0], 'k0p3_psi': point(0.3)[1],
    }

    # Least attainable xi at fixed footrule, against earlier bounds.
    rows = []
    for psi in (-0.05, -0.1, -0.2, -0.3, -0.4, -0.45, -0.49):
        k = k_for_psi(psi)
        x, y = point(k)
        if abs(y-psi) > 1e-12:
            raise ArithmeticError('root finding failed')
        jensen = rockel_jensen(psi)
        if not (psi*psi < jensen < x):
            raise ArithmeticError(f'bound ordering fails at psi={psi}')
        rows.append({'psi': psi, 'k': k, 'min_xi': x, 'rockel_jensen_bound': jensen,
                     'psi_squared': psi*psi, 'ratio_(xi-psi^2)/|psi|^3': (x-psi*psi)/abs(psi)**3})
    out['table'] = rows

    # Observation (not a theorem): (xi - psi^2)/|psi|^3 along the curve.
    ks = np.geomspace(0.05, 2000, 60)
    ratios = []
    for k in ks:
        x, y = point(float(k))
        ratios.append((x-y*y)/abs(y)**3)
    ratios = np.array(ratios)
    out['observation_ratio_xi_minus_psi2_over_abs_psi3'] = {
        'k_range': [float(ks[0]), float(ks[-1])], 'min': float(ratios.min()),
        'max': float(ratios.max()), 'nondecreasing_in_k': bool(np.all(np.diff(ratios) > -1e-9))}

    # Stochastically decreasing block family B_l (Proposition 8).
    fam = []
    for l in (0.0, 0.25, 0.5, 0.75, 1.0):
        xi, psi = block_family(l)
        if abs(xi-(1-l*l)) > 2e-3 or abs(psi+(1-l*l)/2) > 2e-3:
            raise ArithmeticError(f'block family mismatch at l={l}: {xi},{psi}')
        fam.append({'l': l, 'xi_quadrature': xi, 'xi_formula': 1-l*l,
                    'psi_quadrature': psi, 'psi_formula': -(1-l*l)/2})
    out['sd_block_family'] = fam

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2))
    print(json.dumps(out['constants'], indent=2))
    for r in rows:
        print(f"psi={r['psi']:+.2f}  k={r['k']:.6f}  min xi={r['min_xi']:.6f}  "
              f"Jensen={r['rockel_jensen_bound']:.6f}  psi^2={r['psi_squared']:.4f}")
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()

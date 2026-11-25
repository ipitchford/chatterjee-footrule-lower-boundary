"""Independent numerical checks for the (xi,psi) lower rank-copula boundary.

Run with Python 3 and NumPy/SciPy. All reported checkerboard statistics integrate
the cell densities exactly, up to round-off. The certificate is a transportation
linear oracle applied to the gradient of the finite-dimensional convex objective.
"""
import argparse
import csv
import sys

import numpy as np
from scipy.linalg import solve_banded, matmul_toeplitz
from scipy.optimize import isotonic_regression, linear_sum_assignment
from scipy.sparse.linalg import LinearOperator, cg


def clip_primitive(x):
    return np.where(x < 1, x*x/2, x-0.5)


def curve(k, panels):
    n = panels
    t = (np.arange(n)+0.5)/n
    column = np.maximum(1-np.arange(n)/(n*k), 0)/n

    def apply(x):
        return x + matmul_toeplitz(column, x)

    matrix = LinearOperator((n,n), matvec=apply, dtype=float)
    forcing = (clip_primitive(t+k)-clip_primitive(t))/k
    r, info = cg(matrix, forcing, rtol=1e-13, atol=0)
    if info != 0:
        raise RuntimeError(f"linear solve exited with code {info}")
    edges = np.arange(n+1)/n
    integral = np.r_[0., r.cumsum()/n]

    def A(s):
        return np.clip(s,0,1) - np.interp(s,edges,integral) + np.interp(s-k,edges,integral)

    low, high = A(t), A(t+k)
    if min(r.min(), 1-r.max(), np.diff(r).min(), (r-low).min(), (high-r).min()) < -1e-9:
        raise ArithmeticError("branch or monotonicity check failed")
    mean_v_given_t = (1-r)*low + r*high
    xi = 1 - 6*np.mean((2*t-1)*mean_v_given_t)
    # Integrate |u-v| over [0,r] (high branch) and [r,1] (low branch).
    psi = 1 - 3*np.mean(0.5 - r*r + r*high - (1-r)*low)
    return xi, psi


def board(k, n, tol=3e-10):
    m = n-1
    index = np.arange(n)
    dist = abs(index[:,None]-index[None,:]) + np.eye(n)/3
    # h[i,j] is the conditional distribution at interior knot (j+1)/n.
    # After multiplying J by n^2 its Hessian in each row has diagonal 8
    # and off-diagonal 2. The baseline linear term includes endpoint h(1)=1.
    linear = 6*k*np.diff(dist,axis=1)
    linear[:,-1] += 2
    uniform = (1+np.arange(m))/n
    z = np.tile(uniform,(n,1))
    dual = np.zeros_like(z)
    rho = 16
    bands = np.zeros((3,m))
    bands[0,1:] = 2
    bands[1,:] = 8+rho
    bands[2,:-1] = 2
    for it in range(100_000):
        rhs = rho*(z-dual)-linear
        centered = rhs-np.mean(rhs,axis=0)
        h = solve_banded((1,1),bands,centered.T).T + uniform
        old_z = z
        z = np.clip([isotonic_regression(row).x for row in h+dual],0,1)
        dual += h-z
        if max(abs(h-z).max(),rho*abs(z-old_z).max()) < tol:
            break
    else:
        raise RuntimeError("checkerboard solve did not converge")
    full = np.column_stack([np.zeros(n),h,np.ones(n)])
    min_increment = np.diff(full,axis=1).min()
    # h has the correct margins. Mixing with independence restores exact
    # row monotonicity before the convex optimality certificate is evaluated.
    if min_increment < 0:
        mix = -min_increment/(1/n-min_increment)+1e-14
        h = (1-mix)*h+mix*uniform
        full = np.column_stack([np.zeros(n),h,np.ones(n)])
    delta = np.diff(full,axis=1)
    left, right = full[:,:-1], full[:,1:]
    xi = 2*np.mean(left*left+left*right+right*right)-2
    psi = 1-3*np.sum(dist*delta)/(n*n)
    gradient = 8*h + linear
    gradient[:,:-1] += 2*h[:,1:]
    gradient[:,1:] += 2*h[:,:-1]
    # dJ/d(delta) is the reverse cumulative sum of dJ/dh.
    costs = np.column_stack([np.cumsum(gradient[:,::-1],axis=1)[:,::-1],np.zeros(n)])
    rows, cols = linear_sum_assignment(costs)
    gap = (np.sum(costs*delta)-np.sum(costs[rows,cols]))/(n*n)
    if not np.max(abs(delta.sum(axis=0)-1)) < 1e-9:
        raise ArithmeticError('column margins failed')
    if not (np.max(abs(delta.sum(axis=1)-1)) < 1e-9 and delta.min() >= -1e-12):
        raise ArithmeticError('row margins or nonnegativity failed')
    if not -1e-11 <= gap <= 5e-7:
        raise ArithmeticError(f'first-order gap out of range: k={k} n={n} gap={gap}')
    return xi, psi, max(0,gap)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k', nargs='+', type=float, default=[.125,.25,.5,1,2,4,8])
    parser.add_argument('--panels', type=int, default=4096)
    parser.add_argument('--grids', nargs='+', type=int, default=[16,32,64])
    args = parser.parse_args()
    writer = csv.writer(sys.stdout,lineterminator='\n')
    writer.writerow(['k','n','continuum_xi','continuum_psi','continuum_J',
                     'grid_xi','grid_psi','grid_J','grid_minus_continuum_J','grid_FW_gap'])
    for k in args.k:
        if k <= 0:
            parser.error('k must be positive')
        xi, psi = curve(k,args.panels)
        J = xi + 2*k*psi
        for n in args.grids:
            x, p, gap = board(k,n)
            grid_J = x+2*k*p
            writer.writerow([f'{v:.11g}' for v in (k,n,xi,psi,J,x,p,grid_J,grid_J-J,gap)])
            sys.stdout.flush()


if __name__ == '__main__':
    main()

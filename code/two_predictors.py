"""The three-slope example for two independent uniform predictors.

Evaluate the one-dimensional CDF fixed point and compare with a direct
quadratic optimisation on an n x n x n checkerboard. Requires NumPy/SciPy.
"""
import argparse

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import solve_banded
from scipy.optimize import isotonic_regression, linear_sum_assignment


def continuum(k, n=8192):
    if not 0 < k < .25:
        raise ValueError('this check uses the contraction for 0 < k < 1/4')
    t = np.linspace(0,1,n+1)
    s = np.unique(np.r_[t,t+k,t+2*k])
    A = np.clip(s,0,1)
    for _ in range(200):
        H = cumulative_trapezoid(A,s,initial=0)
        h0,h1,h2 = [np.interp(t+j*k,s,H) for j in range(3)]
        a,b = (h1-h0)/k,(h2-h1)/k
        d = b-a
        p0 = (1-a)**2-d*d/2
        p1 = 2*a*(1-b)
        p2 = b*b-d*d/2
        cumul1 = cumulative_trapezoid(1-p0,t,initial=0)
        cumul2 = cumulative_trapezoid(p2,t,initial=0)
        updated = np.clip(s,0,1) - np.interp(s,t,cumul1)+np.interp(s-k,t,cumul1)
        updated += -np.interp(s-k,t,cumul2)+np.interp(s-2*k,t,cumul2)
        err = abs(updated-A).max()
        A = updated
        if err < 1e-13:
            break
    else:
        raise ArithmeticError('fixed point did not converge')
    chosen_sum = a*a*(1-b) + b**3-d*d*(b-d/3)
    slope_Y = np.zeros_like(t)
    xi_integrand = np.zeros_like(t)
    for j,p in enumerate((p0,p1,p2)):
        y = np.interp(t+j*k,s,A)
        slope_Y += p*(t+j*k)*y
        xi_integrand += p*(y-t)**2
    surplus = np.trapezoid(slope_Y-k*chosen_sum,t)
    xi = 6*np.trapezoid(xi_integrand,t)
    objective = 4+4*k-12*surplus
    return xi,(objective-xi)/(4*k),objective


def grid(k,n,tol=5e-10):
    i,j = np.indices((n,n))
    pairs = np.c_[i.ravel(),j.ravel()]
    m = n-1
    ind = np.arange(n)
    D = abs(ind[:,None]-ind)+np.eye(n)/3
    dist = D[pairs[:,0]]+D[pairs[:,1]]
    linear = 6*k*np.diff(dist,axis=1)
    linear[:,-1] += 2
    uniform = (1+np.arange(m))/n
    z = np.tile(uniform,(n*n,1))
    dual = np.zeros_like(z)
    rho = 16
    bands = np.zeros((3,m))
    bands[0,1:] = bands[2,:-1] = 2
    bands[1,:] = 8+rho
    for _ in range(60_000):
        rhs = rho*(z-dual)-linear
        h = solve_banded((1,1),bands,(rhs-rhs.mean(axis=0)).T).T+uniform
        old = z
        z = np.clip([isotonic_regression(row).x for row in h+dual],0,1)
        dual += h-z
        if max(abs(h-z).max(),rho*abs(z-old).max()) < tol:
            break
    else:
        raise ArithmeticError('checkerboard solver failed')
    bp = np.column_stack([np.zeros(n*n),h,np.ones(n*n)])
    min_diff = np.diff(bp,axis=1).min()
    if min_diff < 0:
        w = -min_diff/(1/n-min_diff)+1e-14
        h = (1-w)*h+w*uniform
        bp = np.column_stack([np.zeros(n*n),h,np.ones(n*n)])
    delta = np.diff(bp,axis=1)
    left,right=bp[:,:-1], bp[:,1:]
    xi = 2*np.mean(left**2+left*right+right**2)-2
    pair_psi_sum = 2-3*np.mean(dist*delta)
    objective = xi+2*k*pair_psi_sum
    grad=8*h+linear
    grad[:,:-1] += 2*h[:,1:]
    grad[:,1:] += 2*h[:,:-1]
    cost=np.c_[np.cumsum(grad[:,::-1],axis=1)[:,::-1],np.zeros(n*n)]
    slots=np.repeat(cost,n,axis=1)
    row,col=linear_sum_assignment(slots)
    gap=(np.sum(cost*delta)-np.sum(slots[row,col]))/(n**3)
    if not (max(abs(delta.sum(axis=0)-n))<1e-8 and delta.min()>-1e-12):
        raise ArithmeticError('checkerboard margins or nonnegativity failed')
    if not -1e-11 <= gap <= 5e-7:
        raise ArithmeticError(f'first-order gap out of range: {gap}')
    return xi,pair_psi_sum/2,objective,gap


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k',type=float,default=.125)
    parser.add_argument('--grids',nargs='+',type=int,default=[8,16,24])
    ns=parser.parse_args()
    print('continuum',continuum(ns.k),flush=True)
    for n in ns.grids:
        print('grid',n,grid(ns.k,n),flush=True)

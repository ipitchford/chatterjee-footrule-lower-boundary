"""Finite-interval mode solver for the rank-copula boundary (0 < k <= 1).

The Fredholm kernel has a triangular profile. Twice differentiating reduces
the equation to two finite systems of constant-coefficient ODEs. This program
matches their hyperbolic modes at every breakpoint and imposes the original
Fredholm equation at the two endpoints. Conditional mean and absolute loss
are then integrated by Gauss-Legendre quadrature on each smooth interval.
Requires Python 3, NumPy, and SciPy.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve


def segments(k):
    if not 0 < k <= 1:
        raise ValueError('finite_modes uses 0 < k <= 1')
    m = int(np.floor(1/k+1e-13))
    b = max(1-m*k,0.)
    blocks = [(m+1,0.,b),(m,b,k)] if b > 1e-14 else [(m,0.,k)]
    pieces = []
    offset = 1
    for d, lo, hi in blocks:
        mode = np.arange(1,d+1)
        wave = 2*np.sin(np.pi*mode/(2*(d+1)))/np.sqrt(k)
        for j in range(d):
            eigen = np.sqrt(2/(d+1))*np.sin((j+1)*np.pi*mode/(d+1))
            pieces.append((j*k+lo,j*k+hi,(j+1)/(d+1),offset,wave,eigen))
        offset += 2*d
    return sorted(pieces,key=lambda p:p[0]),offset


def basis(piece, x, operation):
    lo, _, equilibrium, offset, wave, eigen = piece
    z = (x-lo)*wave
    s, c = np.sinh(z), np.cosh(z)
    # cosh(z)-1 is evaluated without cancellation
    cm1 = 2*np.sinh(z/2)**2
    if operation == 'value':
        base, symmetric, antisymmetric = 1.,c,s
    elif operation == 'derivative':
        base, symmetric, antisymmetric = 0.,wave*s,wave*c
    elif operation == 'integral':
        base, symmetric, antisymmetric = x-lo,s/wave,cm1/wave
    elif operation == 'first_moment':
        h = x-lo
        base = h*h/2
        symmetric = h*s/wave-cm1/wave**2
        antisymmetric = h*c/wave-s/wave**2
    else:
        raise ValueError(operation)
    return equilibrium*base, offset, np.r_[eigen*symmetric,eigen*antisymmetric]


def boundary(k, order=24):
    pieces, width = segments(k)

    def vector(piece,x,op='value'):
        constant,start,weights = basis(piece,x,op)
        val = np.zeros(width)
        val[0] = constant
        val[start:start+len(weights)] = weights
        return val

    # Every equation has the form row @ (1,modal_coefficients) = 0.
    eqs = []
    for before,after in zip(pieces,pieces[1:]):
        for op in ('value','derivative'):
            eqs.append(vector(before,before[1],op)-vector(after,after[0],op))
    for t, init, rhs in ((0.,pieces[0],k/2),(1.,pieces[-1],1.)):
        row = vector(init,t)
        row[0] -= rhs
        for p in pieces:
            lo,hi = p[:2]
            if abs((lo+hi)/2-t) < k:
                slope = -1/k if t==0 else 1/k
                start_weight = 1-abs(lo-t)/k
                row += start_weight*vector(p,hi,'integral')+slope*vector(p,hi,'first_moment')
        eqs.append(row)
    E = np.asarray(eqs)
    if E.shape != (width-1,width):
        raise ArithmeticError('mode-matching system has the wrong shape')
    coefficients = np.r_[1.,solve(E[:,1:],-E[:,0])]
    if np.max(abs(E@coefficients)) > 1e-10:
        raise ArithmeticError('mode matching failed')
    right_ends = np.array([p[1] for p in pieces])
    cumulative = np.r_[0.,np.cumsum([vector(p,p[1],'integral')@coefficients for p in pieces])]

    def primitive(s):
        if s <= 0:
            return 0.
        if s >= 1:
            return cumulative[-1]
        j = min(np.searchsorted(right_ends,s),len(pieces)-1)
        return cumulative[j]+vector(pieces[j],s,'integral')@coefficients

    def A(s):
        return np.clip(s,0.,1.)-primitive(s)+primitive(s-k)

    # V lies on two graphs, conditional on the auxiliary T.
    nodes, weights = leggauss(order)
    xi, loss = 0.,0.
    for p in pieces:
        lo,hi = p[:2]
        for node,weight in zip(nodes,weights):
            t, w = (lo+hi+(hi-lo)*node)/2, weight*(hi-lo)/2
            r = vector(p,t)@coefficients
            low, high = A(t),A(t+k)
            xi += 6*w*((1-r)*(low-t)**2+r*(high-t)**2)
            loss += w*(.5-r*r+r*high-(1-r)*low)
    return xi,1-3*loss


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('k',nargs='*',type=float,default=[0.125,0.25,0.3,0.5,0.7,1.])
    for k_value in parser.parse_args().k:
        x,p = boundary(k_value)
        print(f'{k_value:g}, {x:.13f}, {p:.13f}')

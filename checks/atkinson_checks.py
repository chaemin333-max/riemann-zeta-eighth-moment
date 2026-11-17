#!/usr/bin/env python3
"""Consistency checks for the directly proved Gaussian Atkinson bound.

Requires SymPy, NumPy, SciPy and mpmath. The symbolic checks verify algebra;
the quadrature checks catch normalization errors at two finite test points.
The analytic identity and its uniform estimates are proved in the supporting paper, Appendix J.
"""
import json
from fractions import Fraction as Q
from pathlib import Path
import sympy as sp
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma
import mpmath as mp


def symbolic_checks():
    a, x, t = sp.symbols('a x t', positive=True)
    phase = sp.log(1+a/x)/a + x/(1+a)
    assert sp.simplify(sp.diff(phase,x).subs(x,1)) == 0
    assert sp.simplify(sp.diff(phase,x,2).subs(x,1)-(2+a)/(1+a)**2) == 0
    n = t*a*a/(2*sp.pi*(1+a))
    assert sp.simplify(n/(a*sp.diff(n,a))-(1+a)/(2+a)) == 0
    # arsinh(a/(2 sqrt(1+a))) = log(1+a)/2, proved by sinh.
    h=a/(2*sp.sqrt(1+a))
    assert sp.simplify((sp.sqrt(1+a)-1/sp.sqrt(1+a))/2-h) == 0
    root=(a+2)/(2*sp.sqrt(1+a))
    assert sp.simplify(root**2-1-h**2) == 0
    critical=t*sp.log(1+a)+t*a/(1+a)
    atkinson=t*sp.log(1+a)+2*t*h*root-sp.pi*n
    assert sp.simplify(critical-atkinson) == 0
    # Euler remainder has v/12 at the endpoint; Laurent cancellation
    # yields G(u)=1/(2u)-gamma+log(2pi)/2-psi(1-u)/2+sum R.
    u, gam, zeta_prime, psi = sp.symbols('u gam zeta_prime psi')
    F_const=1/(2*u)-gam/2
    beta_zeta_const=zeta_prime+(gam+psi)/2
    G_const=F_const-beta_zeta_const
    assert sp.simplify(G_const-(1/(2*u)-gam-zeta_prime-psi/2)) == 0
    sigmac=Q(27,232)
    exponent=Q(143,116)
    assert 1+2*sigmac==exponent
    assert (94-98*sigmac)/67==exponent
    assert -Q(98,67)<0 and Q(23,19)<exponent
    assert exponent-Q(84599,68640)==Q(509,1990560)
    for sigma in [sigmac,Q(5,38)]:
        assert all(v>=0 for v in [4*(21*sigma-1)/67,(5-38*sigma)/134,
                                  (9-55*sigma)/67,(5-38*sigma)/268])
        assert all(v<0 for v in [(-31+182*sigma)/67,(98-852*sigma)/67,
                                Q(1,8)-Q(9,2)*sigma,1-10*sigma,
                                (36-354*sigma)/67])
        q=(36+48*sigma)/(67*(1-4*sigma))
        assert Q(9,8)<q<=Q(4,3)
    return {'elementary_sixth':str(exponent),'sixth_exponent_cost':str(Q(509,1990560)),
            'stationary_point':True,'curvature':True,
            'parameter_log_derivative':True,'exact_phase':True,
            'Laurent_normalization':True}


def quadrature(T,G,N=256):
    mp.mp.dps=25
    lhs=quad(lambda t:float(abs(mp.zeta(.5+1j*t))**2)*np.exp(-((t-T)/G)**2),
             T-10*G,T+10*G,epsabs=1e-9,limit=400)[0]
    diag=quad(lambda t:(digamma(.5+1j*t).real+2*float(mp.euler)-np.log(2*np.pi))
              *np.exp(-((t-T)/G)**2), T-10*G,T+10*G,epsabs=1e-10,limit=400)[0]
    def real_amplitude(y):
        if y<=0:return 0.
        L=np.log1p(1/y)
        return np.exp(-G**2*L**2/4)*np.cos(T*L)/np.sqrt(y*(1+y))
    divisors=np.zeros(N+1,dtype=int)
    for d in range(1,N+1):divisors[d::d]+=1
    residue=2*np.pi*np.real(np.exp((-.5-1j*T)**2/G**2))
    partial=0.
    rows=[]
    for n in range(1,N+1):
        val,err=quad(real_amplitude,0,np.inf,weight='cos',wvar=2*np.pi*n,
                     epsabs=2e-12,limlst=100,limit=300)
        partial+=int(divisors[n])*val
        if n in (32,64,128,N):
            rhs=diag+4*np.sqrt(np.pi)*G*partial+residue
            rows.append({'N':n,'lhs':lhs,'rhs':float(rhs),'difference':float(lhs-rhs)})
    assert abs(rows[-1]['difference'])<2e-7
    return {'T':T,'G':G,'partials':rows}

if __name__=='__main__':
    result={'symbolic_checks':symbolic_checks(),
            'quadrature_checks':[quadrature(12.,2.),quadrature(30.,3.)],
            'scope':'Algebra and finite numerical consistency; convergence and uniform bounds are proved in the paper.'}
    output=Path(__file__).resolve().parents[1] / 'verification-output' / 'atkinson_checks.json'
    output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result,indent=2))

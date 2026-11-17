"""Exact symbolic identities for sixth-moment/moments.tex and its appendices."""
from fractions import Fraction as Q
import sympy as s

t, u, eta, delta, y = s.symbols('t u eta delta y')
Omega = 2 / (1 + s.sqrt(1 - 3*t))
assert s.simplify(Q(3,4)*t*Omega**2 - Omega + 1) == 0
expansion = s.series(Omega*(4-Omega)/12, t, 0, 4).removeO()
assert s.expand(expansion) == Q(1,4)+t/8+Q(9,64)*t**2+Q(27,128)*t**3

w1, b, c, P = s.symbols('w1 b c P')
w3,w4=w1-b,w1-c
w2=P+w1-b-c
assert s.expand(w1**2+w2**2-w3**2-w4**2
                -(2*b*c+2*P*(w1-b-c)+P**2)) == 0

# The exact external normalization before and after Holder and sampling.
assert Q(3,2)+Q(1,4)==Q(7,4)
assert -Q(1,2)+Q(3,4)+Q(1,4)==Q(1,2)

# Upper-amplitude line from curvature energy: all six affine margins and ranges.
x=s.symbols('sigma');nu=1-4*x;w=(36+48*x)/67;E=(94-500*x)/67
terms=[Q(1,2)+Q(3,2)*nu-w-2*x,
       Q(1,4)+Q(5,4)*nu-w/4-2*x,
       -Q(1,4)+Q(5,4)*nu+w/2-2*x,
       -Q(1,2)+Q(3,2)*nu+Q(3,4)*w-2*x,
       (6+20*nu-4*w-44*x)/17,
       -Q(1,8)+Q(11,8)*nu+w/4-2*x]
margins=[4*(21*x-1)/67,(5-38*x)/134,(9-55*x)/67,0,0,(5-38*x)/268]
for term,margin in zip(terms,margins):
    assert s.simplify(E-term-margin)==0
    for endpoint in (Q(1,8),Q(5,38)):
        assert (E-term).subs(x,endpoint)>=0
conditions=[w-nu/2-Q(1,2),3*nu-1-w,Q(9,8)*nu-1,
            Q(5,2)*nu-Q(3,2),w+nu/2-Q(1,2)-4*x]
for condition in conditions:
    for endpoint in (Q(1,8),Q(5,38)):
        assert condition.subs(x,endpoint)<0
assert s.simplify((E+6*x).subs(x,Q(1,8)))==Q(327,268)<Q(84599,68640)

# Derivative conditions: the three derivative-test lines against
# the j=8 target on all breakpoints. Affinity then gives the full intervals.
k=Q(1,245);ell=Q(269,280)
assert 1+(1-ell)/k==Q(85,8)
for d in range(3,13):
    left=max(Q(2),Q((d-1)**2+1,d))
    right=min(Q(85,8),Q(d*d+1,d+1))
    if left>right: continue
    for tau in (left,right):
        lines=[(tau-d)/Q(d*(d-1)), -Q(1,d*(d-1)),
               -2*tau/Q(d*d*(d-1))]
        target=k*tau+ell-k-1
        assert all(z<=target for z in lines),(d,tau,lines,target)
for tau in (Q(1),Q(2)):
    assert Q(1,2)+tau/6<=k*tau+ell-k
    assert 1-tau/6<=k*tau+ell-k

# Coefficient-box Holder and shift averaging: the H powers must cancel
# at the critical Vinogradov moment. Convert each pair-count term to F.
for d in range(3,13):
    critical=Q(d*(d-1),2)
    h_power=-1+critical-critical*(1-1/critical)-Q(1,2)+Q(1,2)
    assert h_power==0
    D=d*(d-1)
    # (M, lambda M^2, lambda^(-2/d)) after M^(1-2/D),
    # followed by lambda=F M^(-d): pairs are (M exponent, F exponent).
    result=[]
    for m_power,lambda_power in [(1,0),(2,1),(0,-Q(2,d))]:
        result.append((1-Q(2,D)+Q(m_power,D)-d*Q(lambda_power,D),
                       Q(lambda_power,D)))
    assert result==[(1-Q(1,D),0),(1-Q(1,d-1),Q(1,D)),
                    (1,-Q(2,d*d*(d-1)))]

print('PASS: critical-phase identities, duality powers, six upper-range terms,')
print('parameter endpoints, derivative conditions through order 12,')
print('and the shift/mean-value normalization in the derivative proof.')
print('These checks verify symbolic identities and rational inequalities; the paper and cited theorems supply the analytic proofs.')

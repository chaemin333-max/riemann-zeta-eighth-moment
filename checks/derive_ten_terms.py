"""Independent symbolic derivation of the ten terms from the model point-count inequalities.

The ten exponents are generated from point-count inequalities with Z and F_c
left symbolic. SymPy then solves for the intersection. The final assertions compare the result with the paper.
"""
from pathlib import Path
import json
import sympy as S

Q = S.Rational
s, n, w, z, f, r, k, a, d, c = S.symbols('sigma nu w z f r k a d c')
lam = k + (a-k)*n-a*w
mu = c*n+(d-c)*w
K = 1-2*s-w
b = n/2+2*w-Q(3,2)
h = w/2+n/16
ag = 1-Q(3,2)*f-6*s-r/2
A = K+lam/4+2*n+3*w/4-3*z

def transfer(rhs):
    R = S.solve(r-rhs, r)[0]
    before = S.factor(R+z-(1+n)/4)
    # These comparisons use the fixed pairs, so the lower bound for Z may be inserted.
    params = {k:Q(13,414), a:Q(34,207), d:Q(1,56), c:Q(31,280)}
    assert S.diff(before,z).subs(params) < 0
    return S.factor(before.subs(z,f+(1+n)/4))

terms = [
    transfer(K+2*n-2*z),
    transfer(A),
    transfer(K+mu/2+Q(5,2)*n+Q(3,2)*w-5*z),
    transfer(K+mu/2+Q(5,2)*n+Q(3,2)*w-5*z+(n+w-1)/4),
    transfer(A+(k+Q(1,2))*ag/4),
    transfer(A+b/8),
    transfer(K+Q(4,3)*n+w-Q(8,3)*z),
    transfer(K+Q(4,3)*n+w-Q(8,3)*z+(h+(1-n)/2-f-4*s)/3),
    transfer(n+Q(3,4)*w-2*z+Q(3,4)*K+(1+r/2-Q(3,2)*f-6*s)/4),
    transfer(A+b/8+k*ag/4),
]
assert S.simplify(terms[9]+2*f-
    (2*lam+Q(1,2)+2*k+Q(21,2)*n-(16+12*k)*s)/(8+k)) == 0
params = {k:Q(13,414), a:Q(34,207), d:Q(1,56), c:Q(31,280), f:0}
costs = [S.factor(t.subs(params)) for t in terms]
cluster = [S.factor(S.diff(t,f).subs(params)) for t in terms]
assert all(t < 0 for t in cluster)
E4 = (454-2564*s)/279
intersection = S.solve([costs[3].subs(n,1-4*s)-E4,
                   costs[4].subs(n,1-4*s)-E4],(s,w))
sR, wR = intersection[s], intersection[w]
ER = S.factor(E4.subs(s,sR))
p = S.factor(ER+8*sR)
assert p == Q(76024973,51678924)
sL,wL,EL = Q(1,8),Q(5,8),Q(3317,7064)
qL,qR = wL/(1-4*sL),wR/(1-4*sR)
assert sL < sR < Q(5,38) < Q(13,98)
assert Q(175,1448) < sR
assert Q(9,8) < qL < qR < Q(4,3)

endpoints = {}
for name,ss,ww,EE in [('left',sL,wL,EL),('right',sR,wR,ER)]:
    margins = [S.factor(EE-t.subs({s:ss,n:1-4*ss,w:ww})) for t in costs]
    assert all(v >= 0 for v in margins)
    endpoints[name] = dict(zip(map(str,range(1,11)),map(str,margins)))
qq = S.symbols('q')
slopes = [S.factor(S.diff(t.subs(w,qq*n),n)) for t in costs]
assert all(t.subs(qq,qe)>0 for t in slopes for qe in [qL,qR])
assert (ER-EL)/(sR-sL)+8 > 0
assert S.diff(E4+8*s,s) < 0
assert all(p>v for v in [Q(25,17),Q(72,49),EL+8*sL])

# Parameter conditions from eta=N/T and kappa=(TN)^(-1/2).
eta,kap = n-1, -(1+n)/2
bb = eta+kap+2*w
constraints = {
    'W/N^(4/3)':w-Q(4,3)*n,
    'WB/N^(3/2)':w+bb-Q(3,2)*n,
    'B/N^(1/3)':bb-n/3,
    'kappa*W*B':kap+w+bb,
    'B^2/W':2*bb-w,
    'kappa*W':kap+w,
    'eta*N^2/W':eta+2*n-w,
    'eta*N^(1/8)':eta+n/8,
    'Atkinson remainder':Q(5,2)*n-Q(3,2),
    'eta^2*W':2*eta+w,
    'middle pair':w+(n-1)/2-4*s,
}
for key,x in constraints.items():
    for ss,ww in [(sL,wL),(sR,wR)]:
        assert x.subs({s:ss,n:1-4*ss,w:ww}) < 0
    for qe in [qL,qR]:
        slope = S.diff(x.subs(w,qe*n),n)
        assert slope<0 if key=='W/N^(4/3)' else slope>0
absorb = (lam+(k+Q(1,2))*bb-mu-(n+w-1)/2).subs(params)
assert all(absorb.subs({n:nn,w:qe*nn})<0
           for qe in [qL,qR] for nn in [Q(1,3),Q(1,2)])

scalar_curv = transfer(A+(k+Q(1,2))*b/4).subs(params)
scalar_margin=S.factor((E4-scalar_curv.subs(n,1-4*s)).subs(intersection))
assert scalar_margin==Q(2605,432223728)>0
ordinary_curv=transfer(A+b/4).subs(params)
ordinary_excess=S.factor((ordinary_curv.subs(n,1-4*s)-E4).subs(intersection))
assert ordinary_excess==Q(1756631,7780027104)>0

# Direct mixed-term bound: carry the same multiplicities and Z powers.
mixed_prefactor = K+lam/3+Q(13,6)*n+w-Q(11,3)*z+k*b/3
mixed = [transfer(mixed_prefactor),transfer(mixed_prefactor+ag/6)]
assert S.simplify(mixed[0]+Q(8,3)*f-
                  (lam/3+Q(1,12)+Q(5,4)*n-2*s+k*b/3))==0
assert S.simplify(mixed[1]+Q(34,13)*f-
                  (4*lam+3+15*n-36*s+4*k*b)/13)==0
at = {s:sR,w:wR,n:1-4*sR}
mixed_margins = [S.factor(ER-t.subs(params).subs(at)) for t in mixed]
assert all(x>0 for x in mixed_margins)
gap = S.factor((b-(1-6*s-ER/2)).subs(at))
assert gap<0
direct_gain = S.factor((1-params[k])*max(S.Integer(0),gap)/2)
assert direct_gain==0
assert S.diff(costs[3],w)==Q(197,280)>0
assert S.diff(costs[4],w)==-Q(241,883)<0

result = {
    'status':'PASS',
    'scope':'symbolic derivation of the ten terms and exact comparisons',
    'exponent':str(p), 'sigma_R':str(sR), 'w_R':str(wR), 'E_R':str(ER),
    'scalar_curvature_margin':str(scalar_margin),
    'ordinary_curvature_excess':str(ordinary_excess),
    'pair_profile_additional_gain':'0',
    'cluster_exponents':list(map(str,cluster)),
    'endpoint_margins':endpoints,
    'direct_mixed_term_margins':list(map(str,mixed_margins)),
    'B_minus_Ag_exponent':str(gap),
    'direct_H_sum_saving_at_intersection':str(direct_gain),
    'additional_final_moment_saving':'0',
}
Path(__file__).with_name('derived_ten_term_checks.json').write_text(
    json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))

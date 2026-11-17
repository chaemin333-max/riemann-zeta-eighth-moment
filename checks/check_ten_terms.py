"""Exact arithmetic for paper/moments.tex.

Scope: the ten terms at all lengths and amplitudes, and the parameter constraints.
"""
from fractions import Fraction as Q
from pathlib import Path
import json

kp, ap = Q(13,414), Q(34,207)
kd, ad = Q(1,56), Q(31,280)
theta = kp + Q(1,2)
sL, wL, EL = Q(1,8), Q(5,8), Q(3317,7064)
sR, wR, ER = Q(9041221,68905232), Q(65343017,103357848), Q(21777647,51678924)
p = Q(76024973,51678924)
qL, qR = wL/(1-4*sL), wR/(1-4*sR)
p_comparison = Q(346,235)

def costs(s,n,w):
    lam=kp+(ap-kp)*n-ap*w
    mu=ad*n+(kd-ad)*w
    return [
        Q(1,2)+3*n/2-w-2*s,
        lam/4+Q(1,4)+5*n/4-w/4-2*s,
        mu/2-Q(1,4)+5*n/4+w/2-2*s,
        mu/2-Q(1,2)+3*n/2+3*w/4-2*s,
        (4*lam+6+4*kp+20*n-4*w-(44+24*kp)*s)/(17+2*kp),
        lam/4+Q(1,16)+21*n/16-2*s,
        Q(1,3)+2*n/3-2*s,
        Q(1,2)+25*n/48+w/6-10*s/3,
        4*(1+n-6*s)/7,
        (2*lam+Q(1,2)+2*kp+21*n/2-(16+12*kp)*s)/(8+kp),
    ]

def scalar_curvature(s,n,w):
    return costs(s,n,w)[1]+theta*(n/2+2*w-Q(3,2))/4

def comparisons(s,n,w):
    lam=kp+(ap-kp)*n-ap*w
    mu=ad*n+(kd-ad)*w
    b=n/2+2*w-Q(3,2)
    return {
        'W_over_N_four_thirds':w-4*n/3,
        'WB_over_N_three_halves':w+b-3*n/2,
        'B_over_N_one_third':b-n/3,
        'kappa_W_B':3*w-2,
        'translation_B2_over_W':2*b-w,
        'lambda_Btheta_over_mu_sqrtU':lam+theta*b-mu-(n+w-1)/2,
        'kappa_W':w-(1+n)/2,
        'eta_N2_over_W':3*n-1-w,
        'eta_N_one_eighth':9*n/8-1,
        'Atkinson_higher_phase_error':5*n/2-Q(3,2),
        'eta_squared_W':2*n+w-2,
        'refined_pair_second_middle':w+(n-1)/2-4*s,
    }

def fmt(x):
    return {'exact':str(x),'decimal':float(x)}

assert p==ER+8*sR
assert ER==(454-2564*sR)/279
assert qL==Q(5,4) and qL<qR<Q(4,3)
assert Q(175,1448)<sR<Q(13,98)
assert Q(1,8)<sR
assert p_comparison>p

ends={}
for name,s,w,E in [('left',sL,wL,EL),('right',sR,wR,ER)]:
    n=1-4*s
    margins=[E-v for v in costs(s,n,w)]
    assert all(v>=0 for v in margins)
    cmp=comparisons(s,n,w)
    assert all(v<0 for v in cmp.values())
    ends[name]={
        'sigma':fmt(s),'nu_max':fmt(n),'w':fmt(w),'E':fmt(E),
        'term_margins':{str(i+1):fmt(v) for i,v in enumerate(margins)},
        'strict_constraints':{k:fmt(v) for k,v in cmp.items()},
    }
assert costs(sR,1-4*sR,wR)[3]==ER
assert costs(sR,1-4*sR,wR)[4]==ER
scalar_margin=ER-scalar_curvature(sR,1-4*sR,wR)
assert scalar_margin==Q(2605,432223728)>0
for ss,ww,ee in [(sL,wL,EL),(sR,wR,ER)]:
    assert scalar_curvature(ss,1-4*ss,ww)<ee

slopes={}
for name,q in [('left',qL),('right',qR)]:
    ss=[b-a for a,b in zip(costs(Q(0),Q(0),Q(0)),costs(Q(0),Q(1),q))]
    assert all(v>0 for v in ss)
    slopes[name]=[fmt(v) for v in ss]
    # Other non-absorption costs increase with nu. The W/N^(4/3)
    # cost has negative slope and is strictly negative for every nu>0.
    c0=comparisons(Q(0),Q(0),Q(0))
    c1=comparisons(Q(0),Q(1),q)
    for key in c0:
        if key=='W_over_N_four_thirds':
            assert c1[key]-c0[key]<0
            assert c0[key]==0
        elif key!='lambda_Btheta_over_mu_sqrtU':
            assert c1[key]-c0[key]>0

absorption=[]
for q in [qL,qR]:
    for n in [Q(1,3),Q(1,2)]:
        val=comparisons(Q(0),n,q*n)['lambda_Btheta_over_mu_sqrtU']
        assert val<0
        absorption.append(fmt(val))

small_m_pair=(ap-2)/8
small_m_diagonal=(ad-2-9*kd)/8
assert small_m_pair<0 and small_m_diagonal<0
assert Q(1,3)-4*sL<0

moment_slope=(ER-EL)/(sR-sL)+8
assert moment_slope>0 and -Q(2564,279)+8<0
tails={
    'low_amplitude':p-Q(25,17),
    'left_middle_endpoint':p-(EL+8*sL),
    'twelfth_moment_tail':p-Q(72,49),
}
assert all(v>0 for v in tails.values())

# Ordinary curvature exponent and the scalar thin-energy margin at the intersection.
n=1-4*sR
ordinary_curvature=costs(sR,n,wR)[1]+(n/2+2*wR-Q(3,2))/4
ordinary_excess=ordinary_curvature-ER
assert ordinary_excess==Q(1756631,7780027104)>0
diagonal_w_slope=(kd-ad)/2+Q(3,4)
pair_w_slope=-(4*ap+4)/(17+2*kp)
assert diagonal_w_slope==Q(197,280)>0
assert pair_w_slope==-Q(241,883)<0

# The tenth term has cluster multiplicity exponent -2.
assert ((kp-16)-3*kp)/(8+kp)==-2
assert Q(1,4)-Q(3,16)==Q(1,16)
assert Q(5,4)+Q(1,16)==Q(21,16)
assert -Q(1,4)+Q(1,4)==0

lam=kp+(ap-kp)*n-ap*wR
b=n/2+2*wR-Q(3,2)
ag=1-6*sR-ER/2
assert b-ag<0

# Comparison: the ordinary curvature term in place of the two thin terms e6 and e10.
# At nu=1-4 sigma, take the least value over w of the largest term and find where it
# meets the upper-amplitude line (454-2564 sigma)/279: the feasible vertex with largest sigma.
def ordinary_costs(s,n,w):
    c=costs(s,n,w)
    lam=kp+(ap-kp)*n-ap*w
    return c[:5]+c[6:9]+[lam/4-Q(1,8)+11*n/8+w/4-2*s]

def crossing(f):
    base=f(Q(0),Q(1),Q(0))
    ds=[x-y for x,y in zip(f(Q(1),Q(-3),Q(0)),base)]
    dw=[x-y for x,y in zip(f(Q(0),Q(1),Q(1)),base)]
    cons=[(c0-Q(454,279),cs+Q(2564,279),cw) for c0,cs,cw in zip(base,ds,dw)]
    best=None
    for i in range(len(cons)):
        for j in range(i+1,len(cons)):
            (a0,a1,a2),(b0,b1,b2)=cons[i],cons[j]
            det=a1*b2-a2*b1
            if det==0:
                continue
            s=(-a0*b2+a2*b0)/det
            w=(-a1*b0+a0*b1)/det
            if all(c0+cs*s+cw*w<=0 for c0,cs,cw in cons) and (best is None or s>best[0]):
                best=(s,w)
    return best

assert crossing(costs)==(sR,wR)
s_ord,w_ord=crossing(ordinary_costs)
p_ord=(454-2564*s_ord)/279+8*s_ord
assert (s_ord,w_ord)==(Q(21170825,161443546),Q(50992059,80721773))
assert p_ord==Q(118757448,80721773)
assert p_ord-p==Q(35137889893,379237670182932)>0

result={
    'status':'PASS',
    'scope':'exact arithmetic for the ten terms; the analytic proof is in the paper',
    'new_exponent':fmt(p),
    'comparison_exponent':fmt(p_comparison),
    'gain_over_comparison_exponent':fmt(p_comparison-p),
    'pair_profile_additional_gain':fmt(Q(0)),
    'scalar_curvature_margin':fmt(scalar_margin),
    'q_interval':[fmt(qL),fmt(qR)],
    'endpoints':ends,'nu_slopes':slopes,
    'absorption_rectangle_corners':absorption,
    'small_M_pair_margin':fmt(small_m_pair),
    'small_M_diagonal_margin':fmt(small_m_diagonal),
    'middle_moment_slope':fmt(moment_slope),
    'tail_margins':{k:fmt(v) for k,v in tails.items()},
    'ordinary_curvature_excess_at_intersection':fmt(ordinary_excess),
    'B_over_pair_transition_exponent':fmt(b-ag),
    'ordinary_curvature_sigma':fmt(s_ord),
    'ordinary_curvature_exponent':fmt(p_ord),
    'thin_terms_saving':fmt(p_ord-p),
}
Path(__file__).with_name('ten_term_checks.json').write_text(
    json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:result[k] for k in [
    'status','new_exponent','gain_over_comparison_exponent','scalar_curvature_margin',
    'ordinary_curvature_excess_at_intersection']},indent=2))


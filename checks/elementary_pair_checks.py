"""Exact arithmetic for the elementary-pair eighth-moment proof in supplement/moments.tex.
This verifies rational identities and ranges; the supplement supplies the analytic proofs.
"""
from fractions import Fraction as Q
from pathlib import Path
import json

k, a = Q(1,30), Q(1,6)
theta=k+Q(1,2)
sL,wL,EL=Q(1,8),Q(5,8),Q(481,1024)
sR,wR,ER=Q(31801,242699),Q(153305,242699),Q(102674,242699)
p=Q(357082,242699)
qL,qR=wL/(1-4*sL),wR/(1-4*sR)

def costs(s,n,w):
    lam=k+(a-k)*n-a*w
    mu=a*n+(k-a)*w
    return [Q(1,2)+3*n/2-w-2*s,
        lam/4+Q(1,4)+5*n/4-w/4-2*s,
        mu/2-Q(1,4)+5*n/4+w/2-2*s,
        mu/2-Q(1,2)+3*n/2+3*w/4-2*s,
        (4*lam+6+4*k+20*n-4*w-(44+24*k)*s)/(17+2*k),
        lam/4+Q(1,16)+21*n/16-2*s,
        Q(1,3)+2*n/3-2*s,
        Q(1,2)+25*n/48+w/6-10*s/3,
        4*(1+n-6*s)/7,
        (2*lam+Q(1,2)+2*k+21*n/2-(16+12*k)*s)/(8+k)]

def constraints(s,n,w):
    lam=k+(a-k)*n-a*w; mu=a*n+(k-a)*w
    b=n/2+2*w-Q(3,2)
    return {'W_N4/3':w-4*n/3,'WB_N3/2':w+b-3*n/2,
        'B_N1/3':b-n/3,'kappaWB':3*w-2,'B2_W':2*b-w,
        'absorption':lam+theta*b-mu-(n+w-1)/2,
        'kappaW':w-(1+n)/2,'etaN2_W':3*n-1-w,
        'etaN1/8':9*n/8-1,'Atkinson_remainder':5*n/2-Q(3,2),
        'eta2W':2*n+w-2,'middle_pair':w+(n-1)/2-4*s}

assert p==ER+8*sR
assert ER==(242-1340*sR)/157
assert costs(sR,1-4*sR,wR)[3]==ER==costs(sR,1-4*sR,wR)[4]
assert Q(1,8)<sR<Q(9,68)
assert qL==Q(5,4)<qR==Q(30661,23099)<Q(4,3)
ends={}
for name,s,w,E in [('left',sL,wL,EL),('right',sR,wR,ER)]:
    n=1-4*s
    margins=[E-x for x in costs(s,n,w)]
    cond=constraints(s,n,w)
    assert min(margins)>=0 and max(cond.values())<0
    ends[name]={'sigma':str(s),'w':str(w),'E':str(E),
        'margins':list(map(str,margins)),
        'constraints':{t:str(v) for t,v in cond.items()}}
slopes={}
for q in [qL,qR]:
    ss=[y-x for x,y in zip(costs(Q(0),Q(0),Q(0)),costs(Q(0),Q(1),q))]
    assert min(ss)>0
    slopes[str(q)]=list(map(str,ss))
    c0=constraints(Q(0),Q(0),Q(0)); c1=constraints(Q(0),Q(1),q)
    for name in c0:
        if name=='W_N4/3': assert c0[name]==0 and c1[name]<0
        elif name!='absorption': assert c1[name]-c0[name]>0
absorption=[]
for q in [qL,qR]:
    for n in [Q(1,3),Q(1,2)]:
        x=constraints(Q(0),n,q*n)['absorption']
        assert x<0
        absorption.append(str(x))
assert (a-2)/8<0 and (a-2-9*k)/8<0
assert Q(1,3)-4*sL<0
assert ((k-16)-3*k)/(8+k)==-2
assert (ER-EL)/(sR-sL)+8>0

# Derivative orders 3,4,5 cover all nontrivial F >= M.
derivative_ranges=[(3,Q(1),Q(5,2)),(4,Q(5,2),Q(25,8)),(5,Q(25,8),Q(5))]
for d,lo,hi in derivative_ranges:
    alpha=Q(1,2**d-2); beta=1-Q(1,2**(d-2))
    for tau in [lo,hi]:
        target=Q(5,6)+tau/30
        assert 1+(tau-d)*alpha<=target
        assert beta-(tau-d)*alpha<=target
assert Q(5,6)+Q(5,30)==1
# Optimizing the shift and the out-of-range trivial case in induction.
for d in [3,4,5]:
    ao=Q(1,2**(d-1)-2); bo=1-Q(1,2**(d-3))
    an=ao/(2+2*ao); bn=(1+bo)/2
    assert an==Q(1,2**d-2) and bn==1-Q(1,2**(d-2))
    assert bn+Q(2**(d-1)-1,2**d-2)>=1

# The upper-amplitude line (242-1340 sigma)/157 is derived from sixth-power duality.
def upper_costs(s,n,w):
    return [Q(1,2)+3*n/2-w-2*s,
        Q(1,6)+3*n/2-w/3-2*s,
        -Q(3,4)+7*n/4+w-2*s,
        -1+2*n+5*w/4-2*s,
        (6+36*n-8*w-60*s)/25]
upper={}
for s in [sR,Q(9,68)]:
    n=1-4*s; w=(68+184*s)/157; q=w/n; E=(242-1340*s)/157
    margins=[E-x for x in upper_costs(s,n,w)]
    assert min(margins)>=0 and margins[3]==margins[4]==0
    assert Q(9,8)<q<Q(4,3)
    cond=[n-Q(3,5),n/2+2*w-Q(3,2),w-4*n/3,
          3*n-1-w,w-(1+n)/2,2*n+w-2,w+(n-1)/2-4*s]
    assert max(cond)<0
    ss=[y-x for x,y in zip(upper_costs(Q(0),Q(0),Q(0)),upper_costs(Q(0),Q(1),q))]
    assert min(ss)>0
    upper[str(s)]={'w':str(w),'margins':list(map(str,margins)),
        'constraints':list(map(str,cond)),'nu_slopes':list(map(str,ss))}
assert -Q(1340,157)+8==-Q(84,157)<0
# The cluster powers for the four unexpanded upper-amplitude terms are all nonpositive.
assert max([-Q(1),-Q(7,3),-Q(6),-Q(58,25)])<0

# Lower amplitudes use the line derived from curvature energy.
sc=Q(27,232)
assert 1-4*sc==(94-500*sc)/67
for s in [sc,Q(1,8),Q(5,38)]:
    n=1-4*s; w=(36+48*s)/67; q=w/n
    assert Q(9,8)<q<=Q(4,3)
    margins=[4*(21*s-1)/67,(5-38*s)/134,(9-55*s)/67,(5-38*s)/268]
    assert min(margins)>=0
    assert max([(-31+182*s)/67,(98-852*s)/67,Q(1,8)-9*s/2,1-10*s,(36-354*s)/67])<0
assert (94-500*Q(1,8))/67+1==Q(197,134)
assert 1+4*sc<Q(197,134)
assert Q(94,67)+(8-Q(500,67))*Q(1,8)==Q(197,134)
tails={'below_1/8':p-Q(197,134),'left_endpoint':p-(EL+1),'twelfth_tail':p-Q(25,17)}
assert min(tails.values())>0
# Difference from the main exponent of paper/moments.tex.
p_main=Q(76024973,51678924)
assert p-p_main==Q(2428617641,12542423175876)>0
scalar6=costs(sR,1-4*sR,wR)[1]+theta*((1-4*sR)/2+2*wR-Q(3,2))/4
assert scalar6-ER==-Q(494,3640485)<0
result={'status':'PASS',
    'eighth_exponent':str(p),'decimal':float(p),'main_exponent':str(p_main),
    'difference_from_main_exponent':str(p-p_main),'endpoints':ends,'nu_slopes':slopes,
    'absorption_rectangle_corners':absorption,'upper_amplitude':upper,
    'tail_margins':{t:str(v) for t,v in tails.items()},
    'scalar_thin_energy_margin':str(ER-scalar6),
    'cited_inputs':['Iwaniec sparse fourth moment'],
    'proved_in_the_papers':['Gaussian Atkinson short second moment']}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps({t:result[t] for t in ['status','eighth_exponent','decimal','difference_from_main_exponent','cited_inputs','proved_in_the_papers']},indent=2))

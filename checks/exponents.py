"""Exact rational bookkeeping for the large-value line and the sixth-moment exponent.

This script re-derives every exponent from the nine monomials displayed in
Theorem 1.3 of the supporting paper (sixth-moment/moments.tex, display (1.2)).  It uses only
Python's Fraction type.  It verifies arithmetic and parameter inequalities.
The paper supplies the analytic arguments.

Run:  python3 checks/exponents.py
Exit status 0 means every assertion passed.
"""

from fractions import Fraction as Fr

# ---------------------------------------------------------------------------
# Exponent pair: Heath-Brown, j = 8 member of
#   k = 2/((j-1)^2 (j+2)),   1 - l = (3j-2)/(j(j-1)(j+2)).
# ---------------------------------------------------------------------------
j = 8
k = Fr(2, (j - 1) ** 2 * (j + 2))
one_minus_l = Fr(3 * j - 2, j * (j - 1) * (j + 2))
l = 1 - one_minus_l
a = 1 - l + k
assert (k, l, a) == (Fr(1, 245), Fr(269, 280), Fr(17, 392))

# Hypotheses on the pair used in Theorem 1.3.
assert 0 < k < Fr(1, 2)
assert a < Fr(1, 2) + k                     # second-derivative range, first condition
assert a / 3 < (k + Fr(1, 2)) / 8           # second-derivative range, second condition
assert a < 2                                # small-M comparison

# ---------------------------------------------------------------------------
# Affine forms in (1, sigma, nu, w):  value = c0 + cs*sigma + cn*nu + cw*w.
# N = T^nu, W = T^w, V = T^sigma.
# ---------------------------------------------------------------------------
def form(c0=0, cs=0, cn=0, cw=0):
    return (Fr(c0), Fr(cs), Fr(cn), Fr(cw))

def add(*fs):
    return tuple(sum(f[i] for f in fs) for i in range(4))

def scale(c, f):
    return tuple(Fr(c) * x for x in f)

T_ = form(c0=1)
N_ = form(cn=1)
W_ = form(cw=1)
V_ = form(cs=1)

# log_T of the two length costs (Section 6 of the paper):
#   lambda = |eta|^{-k} (W/N)^{-a},  |eta| ~ N/T   ->  (T/N)^k (W/N)^{-a}
#   mu     = W^k (W/N)^{-a}
L_lam = add(scale(k, T_), scale(-k, N_), scale(-a, W_), scale(a, N_))
L_mu = add(scale(k, W_), scale(-a, W_), scale(a, N_))

def mono(lam=0, mu=0, t=0, n=0, w=0, v=0, root=1):
    """log_T of (lambda^lam mu^mu T^t N^n W^w V^v)^root."""
    f = add(scale(lam, L_lam), scale(mu, L_mu), scale(t, T_), scale(n, N_),
            scale(w, W_), scale(v, V_))
    return scale(root, f)

H = Fr(1, 2)
terms = {
    1: mono(t=H, n=Fr(3, 2), w=-1, v=-2),
    2: mono(lam=Fr(1, 4), t=Fr(1, 4), n=Fr(5, 4), w=Fr(-1, 4), v=-2),
    3: mono(mu=H, t=Fr(-1, 4), n=Fr(5, 4), w=H, v=-2),
    4: mono(mu=H, t=-H, n=Fr(3, 2), w=Fr(3, 4), v=-2),
    5: mono(lam=4, t=6 + 4 * k, n=20, w=-4, v=-44 - 24 * k, root=1 / (17 + 2 * k)),
    6: mono(lam=Fr(1, 4), t=Fr(-1, 8), n=Fr(11, 8), w=Fr(1, 4), v=-2),
    7: mono(t=Fr(1, 3), n=Fr(2, 3), v=-2),
    8: mono(t=H, n=Fr(25, 48), w=Fr(1, 6), v=Fr(-10, 3)),
    9: mono(t=Fr(4, 7), n=Fr(4, 7), v=Fr(-24, 7)),
}

def at_nu_star(f):
    """Substitute nu = 1 - 4 sigma.  Returns (c0, cs, cw)."""
    c0, cs, cn, cw = f
    return (c0 + cn, cs - 4 * cn, cw)

# ---------------------------------------------------------------------------
# Balance terms 4 and 5 at nu = 1 - 4 sigma to determine w = omega(sigma).
# ---------------------------------------------------------------------------
d0, ds, dw = (x - y for x, y in zip(at_nu_star(terms[4]), at_nu_star(terms[5])))
omega0, omegas = -d0 / dw, -ds / dw          # omega(sigma) = omega0 + omegas*sigma
assert (omega0, omegas) == (Fr(5 * 1692769, 15938321), Fr(5 * 2461308, 15938321))

def eval_at(f, sigma):
    """Evaluate a term at nu = 1-4 sigma, w = omega(sigma)."""
    c0, cs, cw = at_nu_star(f)
    return c0 + cs * sigma + cw * (omega0 + omegas * sigma)

# E(sigma) is affine; recover it from two points and compare with the claim.
E0 = eval_at(terms[4], Fr(0))
E1 = eval_at(terms[4], Fr(1)) - E0
assert (E0, E1) == (Fr(54928, 38969), Fr(-293156, 38969))
E = lambda s: E0 + E1 * s
assert all(eval_at(terms[5], s) == E(s) for s in (Fr(0), Fr(1)))

# Intersection with the fourth-moment line 1 - 4 sigma.
sigma0 = (1 - E0) / (E1 + 4)
assert sigma0 == Fr(15959, 137280)
assert E(sigma0) == 1 - 4 * sigma0
exponent6 = 1 + 2 * sigma0
assert exponent6 == Fr(84599, 68640)

sig_lo, sig_hi = sigma0, Fr(1, 8)

# ---------------------------------------------------------------------------
# (i) Every term is <= E(sigma) on [sigma0, 1/8] at nu = nu_star.
#     E - e_j is affine in sigma, so the two endpoints suffice.
# ---------------------------------------------------------------------------
margins = {}
for idx, f in terms.items():
    m_lo, m_hi = E(sig_lo) - eval_at(f, sig_lo), E(sig_hi) - eval_at(f, sig_hi)
    assert m_lo >= 0 and m_hi >= 0, idx
    margins[idx] = (m_lo, m_hi)
assert margins[4] == (0, 0) and margins[5] == (0, 0)
for idx in (1, 2, 3, 6, 7, 8, 9):
    assert min(margins[idx]) > 0, idx

# ---------------------------------------------------------------------------
# (ii) W = N^q with q = omega(sigma)/nu_star fixed for all N.  Every term must be
#      increasing in nu for fixed q, so that nu = nu_star is the worst case.
#      d/dnu of a term = cn + q*cw, linear in q: check at q_min and q_max.
# ---------------------------------------------------------------------------
def q_of(s):
    return (omega0 + omegas * s) / (1 - 4 * s)

q_lo, q_hi = q_of(sig_lo), q_of(sig_hi)
# q is increasing in sigma since omegas + 4*omega0 > 0.
assert omegas + 4 * omega0 > 0
assert (q_lo, q_hi) == (Fr(8714085, 7509649), Fr(20004325, 15938321))
assert Fr(9, 8) < q_lo < q_hi < Fr(4, 3)
for idx, f in terms.items():
    _, _, cn, cw = f
    for q in (q_lo, q_hi):
        assert cn + q * cw > 0, (idx, q)

# ---------------------------------------------------------------------------
# (iii) Parameter conditions (Sections 7 and 8 of the paper).  Each quantity is
#       T^{(affine in nu)} with positive nu-slope for fixed q, so nu = nu_star is
#       worst; at nu_star each is affine in sigma, so endpoints suffice.
#       Quantities (log_T), with w = q nu:
#         |kappa| W        : w - nu/2 - 1/2
#         |eta| N^2 / W    : 3 nu - w - 1
#         |eta| N^{1/8}    : (9/8) nu - 1
#         Atkinson remainder N^{5/2} T^{-3/2} : (5/2) nu - 3/2
#         eta^2 W          : 2 nu + w - 2
#         pair absorption 2: (w + nu/2 - 1/2) - 4 sigma   [W sqrt(N/T) < V^4]
# ---------------------------------------------------------------------------
conds = {
    "|kappa|W": form(c0=-H, cn=-H, cw=1),
    "|eta|N^2/W": form(c0=-1, cn=3, cw=-1),
    "|eta|N^(1/8)": form(c0=-1, cn=Fr(9, 8)),
    "Atkinson remainder": form(c0=Fr(-3, 2), cn=Fr(5, 2)),
    "eta^2 W": form(c0=-2, cn=2, cw=1),
    "pair absorption (W sqrt(N/T) < V^4)": form(c0=-H, cs=-4, cn=H, cw=1),
}
cond_values = {}
for name, f in conds.items():
    _, _, cn, cw = f
    for q in (q_lo, q_hi):
        assert cn + q * cw > 0, name          # increasing in nu
    vals = tuple(eval_at(f, s) for s in (sig_lo, sig_hi))
    assert all(v < 0 for v in vals), (name, vals)
    cond_values[name] = vals
# Pair absorption 1: sqrt(T/N) < V^4, worst at N = T^{1/3}: 1/3 - 4 sigma < 0.
assert Fr(1, 3) - 4 * sig_lo < 0
# Range of N is non-empty: 1/3 < nu_star.
assert Fr(1, 3) < 1 - 4 * sig_hi

# ---------------------------------------------------------------------------
# (iv) Small-M comparison: relative price N^{1/8} (W/N)^{-3} is below lambda and
#      mu once W/N >= N^{1/8} (q >= 9/8).
# ---------------------------------------------------------------------------
assert (a - 2) / 8 < 0
assert (a - 2) / 8 - 9 * k / 8 < 0

# ---------------------------------------------------------------------------
# (v) Assembly of the sixth moment over all amplitudes.
# ---------------------------------------------------------------------------
# low amplitudes, fourth moment: layer exponent 1 + 2 sigma <= 1 + 2 sigma0.
# [sigma0, 1/8]: layer E(sigma) + 6 sigma, decreasing.
assert 6 + E1 < 0
assert E(sigma0) + 6 * sigma0 == exponent6
# [1/8, 5/38]: curvature-only line, rederived in Appendix F.
curved = lambda s: Fr(94, 67) - Fr(500, 67) * s
assert 6 - Fr(500, 67) < 0
assert curved(Fr(1, 8)) + Fr(6, 8) == Fr(327, 268) < exponent6
# [5/38, oo): twelfth moment, N(T,V) << T^2 V^{-12}: layer 2 - 6 sigma.
assert 2 - 6 * Fr(5, 38) == Fr(23, 19) < exponent6

# ---------------------------------------------------------------------------
# (vi) Comparisons.
# ---------------------------------------------------------------------------
comparisons = {
    "Heath-Brown (classical) 5/4": Fr(5, 4) - exponent6,
    "de Faveri-Pandey 37/30": Fr(37, 30) - exponent6,
    "rStation paper 38 v3 12867/10438": Fr(12867, 10438) - exponent6,
}
assert comparisons["de Faveri-Pandey 37/30"] == Fr(19, 22880)
assert comparisons["rStation paper 38 v3 12867/10438"] == Fr(73259, 358232160)

# The eighth-moment layer E + 8 sigma increases on this line.
assert 8 + E1 == Fr(18596, 38969) > 0

if __name__ == "__main__":
    print(f"exponent pair (k, l) = ({k}, {l}),  a = 1 - l + k = {a}")
    print(f"omega(sigma) = {omega0} + {omegas} sigma")
    print(f"E(sigma) = ({E0.numerator * 38969 // E0.denominator} "
          f"- {-E1 * 38969} sigma)/38969")
    print(f"sigma0 = {sigma0}")
    print(f"sixth-moment exponent 1 + 2 sigma0 = {exponent6} = {float(exponent6):.15f}")
    print(f"q range: [{q_lo}, {q_hi}] = [{float(q_lo):.6f}, {float(q_hi):.6f}]")
    print("margins E - e_j at (sigma0, 1/8):")
    for idx, (m0, m1) in margins.items():
        print(f"  term {idx}: {m0}, {m1}")
    print("parameter conditions, log_T value at (sigma0, 1/8):")
    for name, (v0, v1) in cond_values.items():
        print(f"  {name}: {v0}, {v1}")
    print("comparisons (previous exponent minus new exponent):")
    for name, d in comparisons.items():
        print(f"  {name}: {d} = {float(d):.12f}")
    print("d/dsigma (E + 8 sigma) =", 8 + E1, "(> 0: eighth-moment assembly uses the separate amplitude estimates)")
    print("ALL CHECKS PASSED (exact rational arithmetic; the analytic proofs are in the paper)")


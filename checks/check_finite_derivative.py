"""Exact interval comparisons for the stated finite derivative theorem."""
from fractions import Fraction as Q
import json
from pathlib import Path

intervals = [
    (3, Q(2), Q(5, 2)),
    (4, Q(5, 2), Q(17, 5)),
    (5, Q(17, 5), Q(13, 3)),
    (6, Q(13, 3), Q(37, 7)),
    (7, Q(37, 7), Q(31, 5)),
]

def target(tau):
    return 1 + tau / 56 - Q(31, 280)

rows = []
for degree, left, right in intervals:
    assert right < degree
    for tau in (left, right):
        costs = [
            1 + (tau - degree) / (degree * (degree - 1)),
            1 - Q(1, degree * (degree - 1)),
            1 - 2 * tau / (degree * degree * (degree - 1)),
        ]
        margins = [target(tau) - c for c in costs]
        assert min(margins) >= 0, (degree, tau, margins)
        rows.append({
            "derivative_order": degree,
            "tau": str(tau),
            "target_exponent": str(target(tau)),
            "three_costs": [str(c) for c in costs],
            "margins": [str(m) for m in margins],
        })

low = []
for tau in (Q(1), Q(2)):
    classical = Q(1, 2) + tau / 6
    margin = target(tau) - classical
    assert margin >= 0
    low.append({"tau": str(tau), "classical_cost": str(classical),
                "margin": str(margin)})
assert target(Q(31, 5)) == 1
assert all(intervals[i][2] == intervals[i + 1][1]
           for i in range(len(intervals) - 1))

result = {
    "status": "PASS",
    "scope": "Affine endpoint comparisons for the stated analytic derivative input",
    "pair": ["1/56", "127/140"],
    "maximum_derivative_order": 7,
    "classical_interval": ["1", "2"],
    "classical_checks": low,
    "derivative_intervals": [[d, str(a), str(b)] for d, a, b in intervals],
    "endpoint_checks": rows,
    "trivial_length_threshold": "31/5",
    "interval_argument": "Each of the three cost differences is affine in tau",
    "analytic_input": "Heath-Brown Theorem 1; rederived in the supporting paper, Appendix H",
}
out = Path(__file__).with_name("finite_derivative_checks.json")
out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
print("PASS: 30 derivative endpoint inequalities, two classical endpoints, and the length threshold")

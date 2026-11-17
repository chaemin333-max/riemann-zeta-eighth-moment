# Computational checks

17 November 2026.

The checks ran with Python 3.14.6. The symbolic calculations use SymPy 1.14.

| File | Check |
|---|---|
| check_manifest.py | SHA-256 hashes of the release files |
| check_ten_terms.py | Ten terms, endpoint bounds, length ranges and amplitude ranges |
| derive_ten_terms.py | Symbolic derivation of the ten terms and their intersection |
| check_finite_derivative.py | Derivative-test bounds at each endpoint |
| exponents.py | Sixth-moment exponents and parameter ranges |
| identity_checks.py | Phase identities and normalizations |
| elementary_pair_checks.py | Eighth-moment proof with direct derivative estimates |
| atkinson_checks.py | Phase formulas and numerical tests of the Gaussian identity |

All seven mathematical checks passed. The rational outputs use exact fractions.
The Atkinson test compares the Gaussian identity with numerical integration at
two heights. Both errors satisfy the 2e-7 tolerance. Fresh numerical results go
to `verification-output/atkinson_checks.json`.

The papers state the analytic inputs and cite their source theorems.
The proofs include the four chord-count terms, both Poisson signs, both Sobolev
classes, the block translations, and the endpoint and tail bounds. The zeta
reduction keeps the cluster multiplicities and covers all amplitudes.

## Derivations and PDF build

The supporting paper keeps the complete derivations of the parts used from
Paper 38: the curvature energy, exact phase, Morse-Fresnel formula and signed
blocks. The thin-energy and support arguments appear once, in the eighth-moment
paper. The supplement cites those results by number.

Tectonic 0.17.0 built all four files. The eighth-moment paper has 11 pages,
the supporting paper has 30, the supplement has 12, and the summary has 2.
The builds have zero unresolved references and zero overfull-box warnings.
The numbered references between the files are read from the auxiliary files
of the earlier builds.

The release archive was extracted into a fresh directory. Building the four
files in the documented order reproduced the same page counts.

## Files

The exact result files use UTF-8 and LF. `SHA256SUMS` lists the release hashes.
The archive was extracted into a fresh directory and checked with
`python checks/check_manifest.py`.

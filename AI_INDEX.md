# AI index: The least functional dependence at negative footrule (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is
a plain introduction with a roadmap). DOI 10.5281/zenodo.24538662 (concept 10.5281/zenodo.24538661).
Status: unrefereed candidate; internal AI-assisted review only.

## Normalisation

(U, V) uniform ranks of continuous (X, Y); copula C; h_C(u, v) = P(V <= v | U = u) (jointly measurable, right-continuous
and nondecreasing in v). xi(C) = 6 int int (h_C - v)^2 du dv (Chatterjee's coefficient, population version, directional
from U to V; equals Rockel's 6 int int (d_1 C)^2 - 2). psi(C) = 1 - 3 E|U - V| = 6 int_0^1 C(t,t) dt - 2 (Spearman's
footrule, range [-1/2, 1]). rho = 12 E[UV] - 3. Pi = independence, M = min, W = (u + v - 1)_+. For k > 0: the
triangular operator (K_k r)(t) = int_0^1 (1 - |z - t|/k)_+ r(z) dz, b_k(t) = k^-1 int_t^{t+k} min(1, max(0, s)) ds,
threshold r_k = (I + K_k)^-1 b_k, level function A_k, extremiser C_k = Law(U, A_k(T + k 1{U < r_k(T)})) with U, T
independent uniforms; (x_k, y_k) = (xi(C_k), psi(C_k)); m(y) = least xi at footrule y. Multiplier mu = 2k.

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A (main) | For each k > 0, C_k is a copula and the unique minimiser of xi + 2k psi over all copulas; (x_k, y_k) given by explicit integrals; y_k decreases strictly from 0 to -1/2 and x_k increases strictly from 0 to 1/2; the full region is -1/2 <= y <= 1, m(y) <= x <= 1, with m(y) = x_k (y_k = y) on (-1/2, 0); C_k is stochastically decreasing; stability xi + 2k psi - x_k - 2k y_k >= 6 ||h_C - h_{C_k}||_2^2 | Theorem 1; Sections 3-5 | written proof; new: the open arc 0 < x < 1/2, -1/2 < y < 0 (Rockel, JCAM 485 (2026) 117466, eq. (16), left it open); upper edge psi = sqrt(xi), corner (1/2, -1/2), segment at y = -1/2, convexity, closedness, SD attainment of the lower edge and \|psi\| <= sqrt(xi) are Rockel's |
| B | Dual (supporting-line) inequality certifying the relaxed transport problem | Lemma 2; Section 3 | written proof |
| C | Closed form for k >= 1 (psi <= -0.4565437...): r_k(t) = 1/2 + d sinh(d(t - 1/2))/(4D), d = sqrt(2/k), D = cosh(d/2); finite piecewise-exponential evaluation for 0 < k < 1 with a nonsingular system; expansions psi = -sqrt(xi) + xi/4 + O(xi^{3/2}) near 0 and psi = -1/2 + (1/2 - xi)^2 + O((1/2 - xi)^3) near 1/2 | Proposition 3; Section 6 | written proof |
| D | xi + 2 psi >= -0.6042071213...; at psi = -0.4565437066... least xi = 0.3088802920...; exact min of xi + psi = -0.1932281... (Rockel's family: -0.193 to three decimals) | Theorem 1 text; Section 1 | values computed from the closed form and finite modes (results/table.json) |
| E | Every attainable pair has a radially symmetric representative; NQD for y < 0, PQD for y > 0 | Corollary 4 | written proof; parallels Orenday Lares-Rockel for (xi, beta) |
| F | PQD or NQD with psi = 0 implies C(u,u) = u^2 and 2 xi + rho <= 1 | Lemma 5 | written proof |
| G | PQD region: {0 < y <= 1, y^2 <= x <= 1} union {(x, 0): 0 <= x < 1/2}; strip copulas C_n with xi = (1 - 1/n)^2/2 | Theorem 6 | written proof; supremum 1/2 not attained |
| H | NQD at psi = 0: xi in [0, 3/4) exactly; Bertino copula B with rho = -1/2, xi = 1/2 | Proposition 7 | written proof; supremum 3/4 not attained |
| I | SD copulas: psi <= 0 (equality only for Pi); xi = 1 only for W; lower edge m(psi) with the same unique extremisers; block family B_l with xi = 1 - l^2, psi = -(1 - l^2)/2; every pair m(y) <= x <= -2y attained by a radially symmetric SD copula | Proposition 8 | written proof |
| J | SD copulas satisfy xi <= -2 psi (would make the SD region exact) | Conjecture 9 | NOT PROVED: conjecture with numerical evidence (sd_probe.py); two proof routes shown insufficient (Section 7) |
| K | Transport with nonuniform margins: unique maximiser of E\|Y - Y'\| + 2k E\|X - Y\|; several displacement terms | Theorems 10 and 11; Section 8 | written proof (existence by Schauder; uniqueness in law) |
| L | (xi - psi^2)/\|psi\|^3 rises from about 1/2 to about 2 along the computed curve | Section 10 | observation only, not a theorem |

## Evidence

| Check | Script | Output | Scope |
|---|---|---|---|
| Closed form vs finite modes at k = 1 (3.3e-16); midpoint Nystrom (4,096 panels) vs exact at k = 1/8 ... 8 (max 5.4e-8); constants; Table 1 by root finding in k; Rockel's Jensen bound; ratio observation; B_l quadrature | [code/boundary_table.py](code/boundary_table.py) | [results/table.json](results/table.json) | double precision; raises on failure |
| Finite-mode evaluation (Gauss-Legendre order 24 per piece) | [code/finite_modes.py](code/finite_modes.py) | [results/finite-modes-check.csv](results/finite-modes-check.csv) | double precision |
| Checkerboard minimisation of xi + 2k psi at n = 16, 32, 64, seven k; grid values above continuum; excess at k = 1: 0.002044, 0.000510, 0.000125; first-order gaps below 4e-10 | [code/verify_boundary.py](code/verify_boundary.py) | [results/verification.csv](results/verification.csv) | consistency check (grid class is a subclass) |
| Two predictors at k = 1/8: scalar -0.029117 vs grids -0.027072, -0.028497, -0.028833 | [code/two_predictors.py](code/two_predictors.py) | [results/two-predictors.txt](results/two-predictors.txt) | consistency check |
| SD vertex search, grids 10, 16, 24, lambda in {1.5, 1.9, 2, 2.1, 3}, independence plus 40 seeded random starts: no positive value for lambda >= 2 (max 1.3e-15) | [code/sd_probe.py](code/sd_probe.py) | [results/sd_probe.json](results/sd_probe.json) | local search: evidence, not proof |
| Figure 1 | [code/plot_boundary.py](code/plot_boundary.py) | [figures/](figures/) | illustration |

All checks raise an error on failure (no assert), so they remain active under python -O; run_all.sh runs them under
-O. The computations are not part of any proof.

## Limits

Population coefficients only. Closed form only for psi <= -0.4565437...; elsewhere double-precision values cross-checked
by independent methods, not interval-certified. Rockel's earlier family was already numerically close; the advance is
exactness, uniqueness, certificate and stability. Conjecture 9 is not proved. The inequality xi >= psi^2 + |psi|^3/2
is suggested by the ratio observation, not proved. Transport fixed points: existence by compactness; constructive
iteration only under the stated contraction conditions.

## Replay and records

Replay command and cross-platform policy (1e-9 absolute for table.json, verification.csv, finite-modes-check.csv;
1e-6 for two-predictors.txt; sd_probe criterion max_over_lam_ge_2 <= 1e-12): [README.md](README.md),
[ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI
[.github/workflows/replay.yml](.github/workflows/replay.yml) with [.github/ci/compare_output.py](.github/ci/compare_output.py).
Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json).
Bibliography: [refs.bib](refs.bib). Licences: [LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes:
[MANIFEST.sha256](MANIFEST.sha256).

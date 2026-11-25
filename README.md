# The least functional dependence at negative footrule

*The lower boundary of the region of Chatterjee's rank correlation and Spearman's footrule: for every negative
footrule, the least possible Chatterjee coefficient, its unique extremal copula and a certificate, completing Rockel's
region*

Version 0.1.0-candidate. DOI [10.5281/zenodo.24538662](https://doi.org/10.5281/zenodo.24538662) (all versions:
[10.5281/zenodo.24538661](https://doi.org/10.5281/zenodo.24538661)). Release page:
https://evidencepress.org/releases/chatterjee-footrule-lower-boundary/

**Status: unrefereed research candidate.** Author: Anonymous. Text CC0-1.0; code MIT (see [LICENSE.md](LICENSE.md)).

## What this is about

**The question.** When analysts look at two measured quantities, X and Y, they often report more than one summary of
how the two depend on each other. Two summaries are compared here. Both use only the ranks of the data.

- *Chatterjee's coefficient* ξ (Chatterjee 2021, now widely used; its population version is due to Dette, Siburg and
  Stoimenov) asks "how nearly is Y a function of X?". It is 0 under independence and 1 when Y is exactly a function of
  X, whatever the shape of that function.
- *Spearman's footrule* ψ asks "how close are the two rankings?". It is 1 for identical rankings, 0 under
  independence, and goes down to -1/2 when the rankings are as far apart as possible.

Which combinations of the two numbers can actually occur together? This set of pairs is the *exact region* of
(ξ, ψ). Knowing it turns into sharp inequalities between the two summaries: for instance, a reported pair outside the
region is impossible.

**What was known.** Rockel (Journal of Computational and Applied Mathematics 485, 2026) found most of the region: its
upper edge ψ = √ξ, the corner at footrule -1/2 (where the least ξ is 1/2, reached only by a checkerboard copula), the
region for "stochastically increasing" data, and the inequality |ψ| ≤ √ξ for every copula. The curved lower edge was
left open: for a footrule between -1/2 and 0, how small can ξ be? He gave a lower bound and a two-parameter family of
examples that was already numerically close to the answer.

**What is new** (all for population coefficients of continuous variables).

1. **The exact lower edge** (Theorem 1). For every footrule strictly between -1/2 and 0, the paper finds the least
   possible ξ, shows that exactly one copula (joint law of the ranks) attains it, and gives a proof certificate valid
   for all copulas. The extremal law has a simple recipe: the response takes one of two levels, depending on whether
   the predictor's rank is below a threshold, and the threshold solves a linear integral equation with a triangular
   kernel. Copulas that come close to the minimum are also close to the extremiser (a stability estimate).
2. **Concrete values.** For every pair of continuous variables, ξ + 2ψ is never below -0.6042071... . At footrule -0.3
   the coefficient ξ must be at least about 0.110, against 0.097 from the earlier lower bound (Table 1 of the paper).
   The edge has an elementary closed form for footrule at most -0.4565..., a finite calculation elsewhere, and explicit
   expansions at both ends.
3. **Restricted classes** (Section 7). Every attainable pair is attained by a radially symmetric copula. Under positive
   or negative quadrant dependence the full regions are determined; at zero footrule ξ can approach, but never reach,
   1/2 and 3/4 respectively. For stochastically decreasing data the lower edge is the same curve, zero footrule forces
   independence, and every pair between the curve and the line ξ = -2ψ is attained.
4. **Transport extensions** (Section 8) to variables with non-uniform distributions and to objectives with several
   displacement terms.

**Idea of the proof.** Written through conditional quantiles, ξ becomes linear rather than quadratic, so minimising
ξ + 2kψ becomes a transport problem: move a uniform square onto a uniform interval with a payoff made of two straight
pieces. The two pieces give a supporting-line (dual) certificate, and asking the output to be uniform gives the integral
equation for the threshold.

**Limits.**
- Population coefficients only, not sample estimates; applying the bounds to data also requires accounting for
  sampling error.
- Closed form only for footrule at most -0.4565...; elsewhere the boundary is a finite numerical calculation, and the
  table values are double-precision computations cross-checked by independent methods, not interval-certified bounds.
- Rockel's earlier family was already numerically close (for ξ + ψ it gives -0.193 to three decimals; the exact
  minimum is -0.1932281...). The advance is exactness, uniqueness, the proof of optimality and stability.
- For stochastically decreasing data, that ξ ≤ -2ψ is the upper edge is a **conjecture** (Conjecture 9), supported by
  numerical search, not a theorem.

**Credit.** The upper edge, the corner at -1/2, the stochastically increasing region, convexity and closedness of the
region, attainment of the lower edge by a stochastically decreasing copula and the inequality |ψ| ≤ √ξ are Rockel's.
The radially symmetric and quadrant-dependence results parallel those of Orenday Lares and Rockel for Chatterjee's ξ
and Blomqvist's β. Ansari and Rockel used an optimal-transport dual construction for Spearman's rho and footrule. The
Bertino copula is from Fredricks and Nelsen.

## Contents

| File | What it is |
|---|---|
| [manuscript.pdf](manuscript.pdf), [manuscript.tex](manuscript.tex) | The paper (17 pages; pdflatex three times; reads `figures/boundary-figure.pdf`; no BibTeX). Section 1 is a plain introduction with a roadmap. |
| [figures/](figures/) | Figure 1 (region and grid-convergence panels), PDF and PNG, made by `code/plot_boundary.py`. |
| [code/verify_boundary.py](code/verify_boundary.py) | Midpoint Nyström solve of the threshold equation (4,096 panels) and independent minimisation of ξ + 2kψ over n × n checkerboard copulas with an assignment-based optimality gap. Output: `results/verification.csv`. |
| [code/finite_modes.py](code/finite_modes.py) | The finite-mode evaluation of Proposition 3 for k ≤ 1. Output: `results/finite-modes-check.csv`. |
| [code/two_predictors.py](code/two_predictors.py) | The two-predictor example of Section 8 against direct grid optimisation. Output: `results/two-predictors.txt`. |
| [code/boundary_table.py](code/boundary_table.py) | Closed form, finite modes and midpoint solver cross-checked; the constants quoted in the text; Table 1; the ratio observation of Section 10; the block family B_l. Output: `results/table.json`. |
| [code/sd_probe.py](code/sd_probe.py) | Numerical evidence for Conjecture 9: local vertex search over stochastically decreasing checkerboard copulas. Output: `results/sd_probe.json`. |
| [code/plot_boundary.py](code/plot_boundary.py), [code/run_all.sh](code/run_all.sh) | The figure, and the driver that regenerates every result and the figure. |
| [requirements.txt](requirements.txt) | Pinned packages: NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 (figure only). |
| [refs.bib](refs.bib) | Bibliography (9 entries, from Crossref and arXiv records). |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits. |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) | What the internal review found and what was done. |
| [RESEARCH_GATES.json](RESEARCH_GATES.json) | Prior-art map, sources read and limitation triage. |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt) | Release replay and environment. |
| [.github/workflows/replay.yml](.github/workflows/replay.yml), [.github/ci/compare_output.py](.github/ci/compare_output.py) | CI replay on Linux and the comparison policy. |
| [CITATION.cff](CITATION.cff), [LICENSE.md](LICENSE.md), [MANIFEST.sha256](MANIFEST.sha256) | Citation metadata, licences, SHA-256 sums. |

## Replay

From the package root, with Python 3.12 and the pinned packages (about 50 to 70 seconds):

```sh
PYTHON="uv run --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 --with matplotlib==3.10.8 python3" sh code/run_all.sh
# or, in an environment with requirements.txt installed:
sh code/run_all.sh
```

`run_all.sh` runs every check under `python -O`; the scripts never use `assert` and raise an error on any failed check,
so they stay active under optimisation. It overwrites `results/` and `figures/`; run it in a copy and compare:

```sh
python3 -I .github/ci/compare_output.py <copy>/results results
```

**Comparison policy** (declared before tagging; the same in the CI workflow, ENVIRONMENT.txt and REPLAY_RECEIPT.md).
`table.json`, `verification.csv` and `finite-modes-check.csv`: same text outside numbers, every number within 1e-9
absolute. `two-predictors.txt` (iterative grid optimisation; the paper quotes six decimals): every number within 1e-6
absolute. `sd_probe.json` (local vertex search, whose vertex values may differ in the last digits): the run raises no
error and `max_over_lam_ge_2` is at most 1e-12. The regenerated figure is not compared byte for byte. On the release
host (macOS arm64, CPython 3.12.11) all five result files were byte-identical to the recorded ones.

## Assurance

Written proofs; the computations test the formulae and the implementation and are not part of any proof. Internal
AI-assisted review only: a combined five-role Claude review panel, one revision round with a bounded extension, a
verification re-review, a proof audit of the new and rewritten proofs, and clean replays with a negative control. No
Astra/Codex review took place. This is not peer review, independent reproduction or formal verification. See
[REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).

## Cite

Anonymous (2026). *The least functional dependence at negative footrule: the lower boundary of the region of
Chatterjee's rank correlation and Spearman's footrule*. Version 0.1.0-candidate. Evidence Press.
https://doi.org/10.5281/zenodo.24538662

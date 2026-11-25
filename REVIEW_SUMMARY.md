# Review summary (public, 0.1.0-candidate)

This records the internal review of the candidate: who looked, what they raised, and what was done. It is
AI-assisted internal review, not peer review or independent reproduction.

## What ran

1. **Panel (round 1).** Five reviewer roles, editor, methodology (proof audit and recomputation), domain (literature,
   attribution and novelty), perspective (readability for non-specialists) and devil's advocate, run as one combined
   Claude pass. The panel recomputed the curve with an independent Nyström solve (at k = 1, (xi, psi) =
   (0.3088802, -0.4565437); the closed forms agreed to 1e-7), computed xi and psi directly from the conditional
   distribution functions without the conditional-quantile identity, checked the dual gap, the endpoint expansions, the
   strip copulas, the Bertino copula and the purification construction, and replayed the bundled scripts.
2. **One revision round** with a bounded extension: the stochastically decreasing case (Proposition 8 proved;
   the upper edge stated as Conjecture 9 with numerical evidence; two proof routes tried and shown insufficient), a new
   context-first introduction, Table 1, expanded proof steps, and scripts with explicit error exits, pinned packages and
   a driver.
3. **Verification re-review** of the revised manuscript against the first version and the cited sources, and a
   **proof audit** of every new or rewritten proof. Their required items were applied by the lead.
4. **Clean replays** of the bundle and of the candidate by the lead, and of the candidate by the publisher from fresh
   copies, with a negative control.

No second-model (Astra/Codex) adversarial review took place for this release.

## Sources

Rockel's paper was read in both the arXiv preprint and the published version (JCAM 485, 117466, open access), which
still leaves its eq. (16) open. Ansari and Rockel (arXiv:2608.20176 and the related arXiv:2506.15897) and Orenday Lares
and Rockel (arXiv:2606.30033) were read and checked against the statements the paper uses. Chatterjee, Dette et al.,
Chiappori-McCann-Pass, Gozlan et al. and Fredricks-Nelsen are cited for standard definitions or context. arXiv listings
and the citation graph of Rockel's preprint, searched up to 24 November 2026, showed no later solution of the open edge.

## Main objections and dispositions

| Objection | Disposition |
|---|---|
| Correctness | No error found by the panel, the re-review or the proof audit. |
| The published version of Rockel's paper and its inequality \|psi\| <= sqrt(xi) were not credited; no comparison with earlier bounds | The published version is the primary citation; its Corollary 2 is credited, with the one-line Cauchy-Schwarz proof noted in Section 5; Table 1 compares the exact least xi with Rockel's Jensen bound and psi^2. |
| The radial and quadrant results parallel Orenday Lares and Rockel for (xi, beta), uncited | Cited; the results are presented as the analogue, and the new contrast is stated: the zero-footrule sections [0, 1/2) and [0, 3/4) are not attained, whereas the beta zero section is the whole interval. |
| Framing overstated the kind of answer | Abstract and introduction say exactly what it is: a closed form for psi <= -0.4565..., a finite piecewise-elementary evaluation otherwise; Rockel's family was already near-optimal (-0.193 against -0.1932281...); the advance is exact values, certificate, unique extremiser and stability. |
| Hard to read for non-specialists | New introduction (question, what was known, what is new, method, roadmap); notation defined before use; figure and table early; every section opens with its purpose. |
| Three proof steps compressed (relaxed problem to copulas; recovering k from C_k; nonsingularity of the finite-mode system) | Each expanded to a complete argument. |
| Checks vanished under python -O; no pinned environment or driver | All asserts replaced by raised errors; requirements.txt; run_all.sh runs every check under -O; negative control. |
| Re-review: the probe script stopped on a 32 x 32 grid (a solver vertex failed its feasibility check) and the table script stalled on a wide root bracket, so the driver did not finish | The probe now runs on the grids 10, 16 and 24, and the text states exactly those grids; the root bracket is [0.02, 1e4]; the full driver replays cleanly. |
| Re-review: theorem numbering of the published Rockel paper; Rockel's convexity, closedness and SD-attainment result (his Theorem 3.2) not credited; Bertino lower edge of the rho-footrule region | Numbering follows the published version (Remark 2(e), Theorem 3.1, Corollary 2); Theorem 3.2 credited in Section 1; Ansari-Rockel cited for the Bertino lower edge. |
| Proof audit: small gaps (case y = 0 in Proposition 8(4); null sets in the case xi = 1; order of two monotonicity statements; witnesses for the two failed proof routes) | Each closed with one or two sentences; the failed-route numbers (0.975 and about 0.737) were checked independently. |
| Suggested extensions | Stochastically decreasing case: proved except the upper edge (Conjecture 9). The clean inequality xi >= psi^2 + \|psi\|^3/2: not proved, reported only as a numerical observation. Interval-certified tables: left open. |

## Limits that remain

Population coefficients only. Closed form only for psi <= -0.4565...; elsewhere the table values are double-precision
computations cross-checked by independent methods, not interval-certified. Conjecture 9 (the stochastically decreasing
upper edge) is not proved. The novelty check is recorded in RESEARCH_GATES.json.

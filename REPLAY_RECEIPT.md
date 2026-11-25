# Replay receipt (release 0.1.0-candidate, 25 November 2026)

The publisher replayed the approved candidate from fresh copies before adding the DOI line. Only `manuscript.tex`
(the DOI line under the author), the rebuilt `manuscript.pdf` and the release records changed afterwards; the code,
the recorded results, the figure files and `requirements.txt` are byte-identical to the replayed copies (checked by
SHA-256 in the frozen-archive gate).

## Run (fresh copy, cwd = package root, macOS arm64)

`PYTHON="uv run --no-project --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 --with matplotlib==3.10.8 python3" sh code/run_all.sh`

Interpreter printed inside the run: CPython 3.12.11, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 (see
[ENVIRONMENT.txt](ENVIRONMENT.txt)). Exit 0 after 46 s; the last line was `run_all: done`, and boundary_table.py printed
`ALL CHECKS PASSED`.

| Output | Result |
|---|---|
| results/verification.csv | byte-identical to the recorded file |
| results/finite-modes-check.csv | byte-identical |
| results/two-predictors.txt | byte-identical |
| results/table.json | byte-identical |
| results/sd_probe.json | byte-identical (max_over_lam_ge_2 = 1.3e-15) |
| figures/boundary-figure.png | byte-identical |
| figures/boundary-figure.pdf | bytes differ (font subsetting); renders identically at 100 dpi |

The lead's earlier replay of the same candidate (24 November 2026, same pinned packages) also gave byte-identical CSV and
TXT files and JSON-equal table and probe files.

## Negative control (fresh copy, `python3 -I -B -O code/boundary_table.py`)

| Mutation | Result |
|---|---|
| Block family B_l: the expected xi changed from 1 - l^2 to 1 - l^2/2 | `ArithmeticError: block family mismatch at l=0.25`, exit 1 |

## Comparison policy (declared before tagging; also in README.md, ENVIRONMENT.txt and the CI workflow)

- table.json, verification.csv, finite-modes-check.csv: same text outside numbers; numbers within 1e-9 absolute.
- two-predictors.txt: numbers within 1e-6 absolute.
- sd_probe.json: no error and max_over_lam_ge_2 <= 1e-12.

The comparison script is `.github/ci/compare_output.py`; it accepted the release output (all differences 0) and rejected
a copy with one value in verification.csv moved by 2e-9.

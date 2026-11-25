# SPDX-License-Identifier: MIT
"""Compare regenerated result files with the recorded ones under the declared cross-platform policy.

usage: compare_output.py NEW_RESULTS_DIR RECORDED_RESULTS_DIR

Policy (declared before tagging; the same in README.md, ENVIRONMENT.txt and REPLAY_RECEIPT.md):
- table.json, verification.csv, finite-modes-check.csv: same structure and the same text outside numbers; every number
  within 1e-9 absolute.
- two-predictors.txt (iterative grid optimisation; the paper quotes six decimals): same text outside numbers; every
  number within 1e-6 absolute.
- sd_probe.json (local vertex search; vertex values may differ in the last digits): the run raised no error and
  max_over_lam_ge_2 <= 1e-12.
- Interpreter or library version strings are not part of these files and are not compared.
Exits 1 on any difference.
"""
import json
import re
import sys
from pathlib import Path

NUM = re.compile(r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?')


def text_cmp(name, new, old, tol, bad):
    a, b = new.splitlines(), old.splitlines()
    if len(a) != len(b):
        bad.append(f'{name}: line count {len(a)} != {len(b)}')
        return
    worst = 0.0
    for i, (x, y) in enumerate(zip(a, b), 1):
        if NUM.sub('#', x) != NUM.sub('#', y):
            bad.append(f'{name} line {i} differs outside numbers: {x!r} != {y!r}')
            continue
        for p, q in zip(NUM.findall(x), NUM.findall(y)):
            d = abs(float(p) - float(q)); worst = max(worst, d)
            if d > tol:
                bad.append(f'{name} line {i}: {p} vs {q} differs by {d:g} > {tol:g}')
    print(f'{name}: {len(a)} lines, largest numeric difference {worst:g} (tolerance {tol:g})')


def json_cmp(name, a, b, tol, path, bad, acc):
    if isinstance(b, dict):
        if not isinstance(a, dict) or sorted(a) != sorted(b):
            bad.append(f'{name}{path}: keys differ'); return
        for k in b:
            json_cmp(name, a[k], b[k], tol, f'{path}.{k}', bad, acc)
    elif isinstance(b, list):
        if not isinstance(a, list) or len(a) != len(b):
            bad.append(f'{name}{path}: list length differs'); return
        for i, (x, y) in enumerate(zip(a, b)):
            json_cmp(name, x, y, tol, f'{path}[{i}]', bad, acc)
    elif isinstance(b, bool) or b is None or isinstance(b, str):
        if a != b:
            bad.append(f'{name}{path}: {a!r} != {b!r}')
    else:
        d = abs(float(a) - float(b)); acc[0] = max(acc[0], d)
        if d > tol:
            bad.append(f'{name}{path}: {a} vs {b} differs by {d:g} > {tol:g}')


def main():
    new, old = Path(sys.argv[1]), Path(sys.argv[2])
    bad = []
    for f in ('verification.csv', 'finite-modes-check.csv'):
        text_cmp(f, (new / f).read_text(), (old / f).read_text(), 1e-9, bad)
    text_cmp('two-predictors.txt', (new / 'two-predictors.txt').read_text(), (old / 'two-predictors.txt').read_text(), 1e-6, bad)
    acc = [0.0]
    json_cmp('table.json', json.loads((new / 'table.json').read_text()), json.loads((old / 'table.json').read_text()), 1e-9, '', bad, acc)
    print(f'table.json: largest numeric difference {acc[0]:g} (tolerance 1e-09)')
    sd = json.loads((new / 'sd_probe.json').read_text())
    m = sd.get('max_over_lam_ge_2')
    if not isinstance(m, (int, float)) or m > 1e-12:
        bad.append(f'sd_probe.json: max_over_lam_ge_2 = {m!r} (must be <= 1e-12)')
    print(f'sd_probe.json: max_over_lam_ge_2 = {m!r} (criterion <= 1e-12)')
    if bad:
        print('DIFFERENCES:')
        for b in bad:
            print(' ', b)
        sys.exit(1)
    print('OK: all result files agree with the recorded ones under the declared policy')


if __name__ == '__main__':
    main()

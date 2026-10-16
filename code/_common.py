"""Shared helpers: paths, result writing, and pass/fail bookkeeping."""
import json, csv, os, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, 'results')
FIGURES = os.path.join(ROOT, 'figures')
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(FIGURES, exist_ok=True)


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    return o


def write_json(name, data):
    path = os.path.join(RESULTS, name)
    with open(path, 'w') as f:
        json.dump(_clean(data), f, indent=2)
    return path


def write_csv(name, header, rows):
    path = os.path.join(RESULTS, name)
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([_clean(x) for x in r])
    return path


def read_json(name):
    with open(os.path.join(RESULTS, name)) as f:
        return json.load(f)


class Checks:
    """Collects comparisons against the numbers printed in the preprint."""

    def __init__(self, title):
        self.title = title
        self.items = []
        print(f'\n=== {title} ===')

    def close(self, name, value, target, tol, ref=''):
        ok = abs(value - target) <= tol
        self._add(name, ok, f'got {value:.6g}, preprint {target:g} (tol {tol:g}){" - " + ref if ref else ""}')
        return ok

    def less(self, name, value, bound, ref=''):
        ok = value < bound
        self._add(name, ok, f'got {value:.3g} < {bound:g}{" - " + ref if ref else ""}')
        return ok

    def true(self, name, cond, detail=''):
        self._add(name, bool(cond), detail)
        return bool(cond)

    def _add(self, name, ok, detail):
        self.items.append({'check': name, 'pass': ok, 'detail': detail})
        print(f'  [{"PASS" if ok else "FAIL"}] {name}: {detail}')

    def finish(self):
        n_fail = sum(not c['pass'] for c in self.items)
        print(f'--- {self.title}: {len(self.items) - n_fail}/{len(self.items)} checks passed')
        return self.items


def run_main(main):
    """Standard entry point: run, and exit with status 1 if any check fails."""
    t0 = time.time()
    quick = '--quick' in sys.argv
    items = main(quick=quick)
    bad = [c for c in items if not c['pass']]
    print(f'(runtime {time.time() - t0:.1f} s)')
    if bad:
        print(f'\nFAILED: {len(bad)} check(s) did not reproduce the preprint.')
        sys.exit(1)
    print('\nSUCCESS: all checks in this script reproduce the preprint.')

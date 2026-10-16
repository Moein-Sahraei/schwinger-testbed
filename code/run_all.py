"""Run every analysis script, then regenerate the figures.

    python code/run_all.py           full run (reproduces every number checked)
    python code/run_all.py --quick   faster subset used by GitHub Actions

Exits with status 1 if any check fails to reproduce the preprint.
"""
import sys, os, time, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import write_json

SCRIPTS = ['schwinger_ed', 'schwinger_pauli', 'schwinger_trotter', 'schwinger_commutator',
           'schwinger_proper_gap', 'schwinger_circuit', 'schwinger_hierarchy', 'schwinger_ordering',
           'make_figures']


def main():
    quick = '--quick' in sys.argv
    t0 = time.time()
    summary, failed = [], []
    for name in SCRIPTS:
        t = time.time()
        items = importlib.import_module(name).main(quick=quick)
        n_fail = sum(not c['pass'] for c in items)
        summary.append({'script': f'code/{name}.py', 'checks': len(items), 'failed': n_fail,
                        'seconds': round(time.time() - t, 1)})
        failed += [f"{name}: {c['check']} ({c['detail']})" for c in items if not c['pass']]
    total = sum(s['checks'] for s in summary)
    write_json('run_summary.json', {'quick': quick, 'total_checks': total, 'failed': len(failed),
                                    'runtime_seconds': round(time.time() - t0, 1), 'scripts': summary})
    print('\n' + '=' * 72)
    for s in summary:
        print(f"  {s['script']:<32} {s['checks'] - s['failed']:>3}/{s['checks']:<3} checks   {s['seconds']:>6.1f} s")
    print('=' * 72)
    print(f'Total runtime {time.time() - t0:.1f} s, {total - len(failed)}/{total} checks passed'
          + (' (quick mode)' if quick else ''))
    if failed:
        print('\nFAILED - these results did not reproduce the preprint:')
        for f in failed:
            print('  -', f)
        sys.exit(1)
    print('\nSUCCESS: all numbers checked reproduce the preprint, and all figures were regenerated.')


if __name__ == '__main__':
    main()

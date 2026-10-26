"""Batch fitting: one isolated mode and one spectrum per file."""
import argparse
import csv
from pathlib import Path

from qfit import MODELS, fit_resonance, load_spectrum, save_result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    parser.add_argument('--pattern', default='*.txt')
    parser.add_argument('--kind', choices=['amplitude', 'power', 'db', 'complex'], default='amplitude')
    parser.add_argument('--model', choices=['auto', *MODELS], default='fano')
    parser.add_argument('--feature', choices=['auto', 'peak', 'dip'], default='peak')
    parser.add_argument('--window', type=float, nargs=2, required=True)
    parser.add_argument('--frequency-unit', default='THz')
    parser.add_argument('--baseline', choices=['constant', 'linear'], default='linear')
    parser.add_argument('--fcol', type=int, default=0)
    parser.add_argument('--ycol', type=int, default=1)
    parser.add_argument('--imag-col', type=int)
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--output-dir', type=Path, default=Path('results/qfit_batch'))
    args = parser.parse_args(argv)
    if (args.kind == 'complex') != (args.imag_col is not None):
        parser.error('--kind complex and --imag-col must be specified together')
    if args.output_dir.resolve() == args.folder.resolve():
        parser.error('Output directory must differ from the input directory.')
    files = sorted(p for p in args.folder.glob(args.pattern) if p.is_file()
                   and not p.resolve().is_relative_to(args.output_dir.resolve()))
    if not files: parser.error('No files match the input pattern.')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for source in files:
        # Include extension to prevent a.txt / a.csv basename collisions.
        relative = source.relative_to(args.folder)
        prefix = args.output_dir / relative.with_name(source.name + '_fit')
        row = {'file': str(source), 'status': 'failed'}
        try:
            f, y = load_spectrum(source, args.fcol, args.ycol, args.imag_col)
            r = fit_resonance(f, y, kind=args.kind, model=args.model, window=args.window,
                               feature=args.feature, baseline=args.baseline, frequency_unit=args.frequency_unit)
            save_result(r, prefix, plot=not args.no_plots)
            row.update(status='fit_completed', model=r.model, f0=r.f0, linewidth=r.linewidth,
                       q_resonance=r.q_resonance, q_spectral=r.q_spectral, q_3db=r.q_3db,
                       rmse_power=r.diagnostics['rmse_power'], warnings='; '.join(r.warnings))
        except (ValueError, OSError, RuntimeError) as exc:
            row['error'] = str(exc)
        rows.append(row)
        print(row['status'], source.name)
    fields = ['file', 'status', 'model', 'f0', 'linewidth', 'q_resonance', 'q_spectral',
              'q_3db', 'rmse_power', 'warnings', 'error']
    with (args.output_dir / 'summary.csv').open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    # A completed fit with warnings still needs scientific review.
    return 1 if any(r['status'] == 'failed' for r in rows) else 0


if __name__ == '__main__':
    raise SystemExit(main())

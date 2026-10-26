"""Single-resonance power-spectrum fits with explicit linewidth/Q conventions.

Requires numpy, scipy and (only for plots) matplotlib. No signal is smoothed.
Fano uses eps=2*(f-f0)/Gamma, so the fitted pole Q is f0/Gamma.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import re
from typing import Any

import numpy as np
from scipy.optimize import brentq, least_squares
from scipy.special import voigt_profile

MODELS = ("lorentzian", "fano", "gaussian", "voigt", "pseudo_voigt")


def to_power(values, kind="amplitude"):
    """Convert S-parameter magnitude, power, complex amplitude or dB to power."""
    a = np.asarray(values)
    if kind == "complex":
        p = np.abs(a) ** 2
    else:
        if np.iscomplexobj(a):
            raise ValueError("Complex data require kind='complex'.")
        a = a.astype(float)
        if kind in ("amplitude", "magnitude", "power"):
            if np.any(a < 0):
                raise ValueError(f"Negative {kind} values; check the input kind.")
            p = a ** 2 if kind != "power" else a.copy()
        elif kind == "db":
            with np.errstate(over="ignore"):
                p = 10.0 ** (a / 10.0)
        else:
            raise ValueError("kind must be amplitude, power, db or complex")
    if not np.all(np.isfinite(p)):
        raise ValueError("Input contains NaN/Inf or dB conversion overflow.")
    return p


def load_spectrum(path, fcol=0, ycol=1, imag_col=None, delimiter=None, skiprows=0):
    """Read numeric TXT/CSV/TSV, including CST headers. Columns are zero-based.

    Text header/comment lines are ignored. A numeric line with missing columns,
    NaN/Inf, or repeated frequency later fails validation instead of being hidden.
    Decimal separator must be a dot. Concatenated CST sweep blocks are rejected.
    """
    if min(fcol, ycol) < 0 or (imag_col is not None and imag_col < 0):
        raise ValueError("Column indices must be nonnegative.")
    rows = []
    maxcol = max(fcol, ycol, imag_col if imag_col is not None else 0)
    for lineno, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        if lineno <= skiprows:
            continue
        line = line.strip()
        if not line or line.startswith(("#", "!", "//", "%")):
            continue
        parts = re.split(r"[\s,;]+", line) if delimiter is None else line.split(delimiter)
        try:
            float(parts[0].replace("D", "E").replace("d", "e"))
        except (ValueError, IndexError):
            continue
        if len(parts) <= maxcol:
            raise ValueError(f"Numeric row {lineno} has too few columns.")
        try:
            vals = [float(v.replace("D", "E").replace("d", "e")) for v in parts]
        except ValueError as exc:
            raise ValueError(f"Invalid numeric row {lineno}.") from exc
        yy = vals[ycol] if imag_col is None else complex(vals[ycol], vals[imag_col])
        rows.append((vals[fcol], yy))
    if not rows:
        raise ValueError("No numeric data found.")
    f = np.array([r[0] for r in rows], float)
    y = np.array([r[1] for r in rows])
    return f, y


def _profiles(x, shape, model):
    c, w = shape[0], np.exp(shape[1])
    z = (x - c) / w
    lorentz = 1.0 / (1.0 + 4.0 * z ** 2)
    if model == "lorentzian":
        return [lorentz]
    if model == "fano":
        return [lorentz, 2.0 * z * lorentz]
    gaussian = np.exp(-4.0 * np.log(2.0) * z ** 2)
    if model == "gaussian":
        return [gaussian]
    if model == "pseudo_voigt":
        return [(1.0 - shape[2]) * gaussian + shape[2] * lorentz]
    if model == "voigt":
        # w is Gaussian FWHM; exp(shape[2]) is Lorentzian FWHM.
        sigma = w / np.sqrt(8.0 * np.log(2.0))
        gamma = np.exp(shape[2]) / 2.0
        return [voigt_profile(x - c, sigma, gamma) / voigt_profile(0.0, sigma, gamma)]
    raise ValueError(f"Unknown model: {model}")


def _design(x, shape, model, baseline):
    cols = [np.ones_like(x)]
    if baseline == "linear":
        cols.append(x)
    return np.column_stack(cols + _profiles(x, shape, model))


@dataclass
class FitResult:
    model: str
    frequency_unit: str
    f0: float
    linewidth: float
    q_resonance: float | None
    q_spectral: float | None
    feature_frequency: float
    feature_kind: str
    feature_fwhm: float | None
    q_3db: float | None
    bandwidth_3db: float | None
    parameters: dict[str, Any]
    diagnostics: dict[str, Any]
    warnings: list[str] = field(default_factory=list)
    candidates: list[dict[str, Any]] = field(default_factory=list)
    frequency: np.ndarray = field(default_factory=lambda: np.array([]), repr=False)
    power: np.ndarray = field(default_factory=lambda: np.array([]), repr=False)
    fitted_power: np.ndarray = field(default_factory=lambda: np.array([]), repr=False)
    fitted_baseline: np.ndarray = field(default_factory=lambda: np.array([]), repr=False)

    def to_dict(self):
        d = asdict(self)
        for name in ("frequency", "power", "fitted_power", "fitted_baseline"):
            d.pop(name)
        return d


def _nearest_crossings(x, y, index, level):
    """Nearest bracketed crossings on each side; missing sides stay None."""
    def interpolate(a, b):
        if y[b] == y[a]:
            return float((x[a] + x[b]) / 2.0)
        return float(x[a] + (level - y[a]) * (x[b] - x[a]) / (y[b] - y[a]))
    left = right = None
    for a in range(index - 1, -1, -1):
        if (y[a] - level) * (y[a + 1] - level) <= 0:
            left = interpolate(a, a + 1)
            break
    for a in range(index, len(x) - 1):
        if (y[a] - level) * (y[a + 1] - level) <= 0:
            right = interpolate(a, a + 1)
            break
    return left, right


def _one_fit(f, power, model, baseline, feature, frequency_unit, n_starts):
    mid, span = float((f[-1] + f[0]) / 2), float(np.ptp(f))
    x = (f - mid) / span
    offset, scale = float(np.median(power)), float(np.ptp(power))
    y = (power - offset) / scale
    nbase = 2 if baseline == "linear" else 1
    minw, maxw = 1e-6, 2.0
    lo, hi = [-.5, np.log(minw)], [.5, np.log(maxw)]
    if model == "voigt":
        lo.append(np.log(minw)); hi.append(np.log(maxw))
    elif model == "pseudo_voigt":
        lo.append(0.0); hi.append(1.0)

    def project(shape):
        design = _design(x, shape, model, baseline)
        beta = np.linalg.lstsq(design, y, rcond=None)[0]
        return design @ beta - y, beta

    # Variable projection fits the baseline and signed contrast coefficients
    # exactly at each nonlinear center/width step. Both peaks and dips work.
    edge = max(2, len(x) // 10)
    edge_i = np.r_[np.arange(edge), np.arange(len(x) - edge, len(x))]
    b = np.polyfit(x[edge_i], y[edge_i], 1)
    rough = y - np.polyval(b, x)
    center0 = float(x[np.argmax(np.abs(rough))])
    centers = [center0, 0.0, -.25, .25, -.4, .4]
    widths = [.03, .10, .30, .65]
    starts = [(c, w) for w in widths for c in centers]
    if n_starts < 1:
        raise ValueError("n_starts must be positive")
    # Interleave width and center diversity for a small start budget.
    starts = starts[::3] + [s for i, s in enumerate(starts) if i % 3]
    best = None
    for c, w in starts[:n_starts]:
        initial = [np.clip(c, -.499, .499), np.log(w)]
        if model == "voigt": initial.append(np.log(w / 2))
        if model == "pseudo_voigt": initial.append(.5)
        opt = least_squares(lambda s: project(s)[0], initial, bounds=(lo, hi),
                            max_nfev=600, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        rss = float(np.dot(opt.fun, opt.fun))
        if best is None or rss < best[0]: best = (rss, opt)
    _, opt = best
    residual, beta = project(opt.x)
    beta = beta * scale
    beta[0] += offset
    fitted = _design(x, opt.x, model, baseline) @ beta
    base = beta[0] + (beta[1] * x if nbase == 2 else np.zeros_like(x))
    center, width = mid + opt.x[0] * span, np.exp(opt.x[1]) * span
    params = {"baseline_at_window_midpoint": float(beta[0]),
              "baseline_slope_per_frequency_unit": float(beta[1] / span) if nbase == 2 else 0.0,
              "window_midpoint": mid, "input_kind_fitted": "power"}
    q_resonance = None
    if model in ("lorentzian", "fano"):
        q_resonance = float(center / width)
        params["Gamma_FWHM_pole"] = float(width)
    if model == "fano":
        c, d = float(beta[nbase]), float(beta[nbase + 1])
        params.update(lorentz_coefficient=c, dispersive_coefficient=d,
                      epsilon_definition="2*(f-f0)/Gamma")
        norm = np.hypot(c, d)
        if norm == 0 or abs(d) < 1e-8 * norm:
            q = 0.0 if c < 0 else None
        else:
            q = (norm + c) / d if c >= 0 else d / (norm - c)
        params["fano_q_canonical_positive_scale"] = None if q is None else float(q)
    else:
        params["signed_height"] = float(beta[nbase])
    if model == "pseudo_voigt": params["lorentz_fraction"] = float(opt.x[2])
    if model == "voigt":
        wg, wl = width, np.exp(opt.x[2]) * span
        params.update(gaussian_fwhm=float(wg), lorentzian_fwhm=float(wl),
                      q_lorentz_component=float(center / wl))
        sigma, gamma = wg / np.sqrt(8 * np.log(2)), wl / 2
        v0 = voigt_profile(0., sigma, gamma)
        width = 2 * brentq(lambda t: voigt_profile(t, sigma, gamma) / v0 - .5,
                           0., 10 * (wg + wl), xtol=(wg + wl) * 1e-11)

    # Contrast FWHM: subtract fitted local baseline, then use the selected
    # positive/negative lobe. For Fano this differs from the pole linewidth.
    xd = np.linspace(-.5, .5, 40001)
    fd = mid + xd * span
    yd = _design(xd, opt.x, model, baseline) @ beta
    bd = beta[0] + (beta[1] * xd if nbase == 2 else np.zeros_like(xd))
    contrast = yd - bd
    if feature == "auto":
        idx = int(np.argmax(np.abs(contrast)))
        selected = "peak" if contrast[idx] >= 0 else "dip"
    else:
        selected = feature
        idx = int(np.argmax(contrast) if feature == "peak" else np.argmin(contrast))
    strength = float(contrast[idx])
    left = right = None
    if (selected == "peak" and strength > 0) or (selected == "dip" and strength < 0):
        left, right = _nearest_crossings(fd, contrast, idx, strength / 2)
    fwhm = None if left is None or right is None else right - left
    ff = float(fd[idx])
    q_spectral = None if not fwhm or fwhm <= 0 else ff / fwhm
    q3 = bw3 = peak3 = None
    if selected == "peak":
        maxima = np.flatnonzero((yd[1:-1] >= yd[:-2]) & (yd[1:-1] > yd[2:])) + 1
        maxima = maxima[np.abs(fd[maxima] - ff) <= width]
        if len(maxima):
            i3 = int(maxima[np.argmin(np.abs(fd[maxima] - ff))])
            peak3 = float(fd[i3])
            if yd[i3] > 0:
                l3, r3 = _nearest_crossings(fd, yd, i3, yd[i3] / 2)
                if l3 is not None and r3 is not None:
                    bw3 = r3 - l3
                    if bw3 > 0: q3 = peak3 / bw3

    rss = float(np.sum((power - fitted) ** 2))
    k, n = len(opt.x) + len(beta), len(f)
    variance = max(rss / n, np.finfo(float).tiny)
    aic = n * np.log(variance) + 2 * k
    aicc = float(aic + 2 * k * (k + 1) / (n - k - 1)) if n > k + 1 else None
    ss = float(np.sum((power - np.mean(power)) ** 2))
    warnings = []
    if not opt.success: warnings.append("Optimizer did not converge.")
    if model not in ("fano", "lorentzian"):
        warnings.append("Spectral Q is a linewidth descriptor; it is not automatically a cavity decay Q.")
    if model == "fano":
        warnings.append("Fano pole Q=f0/Gamma and feature-width Q have different definitions.")
    if model == "voigt":
        warnings.append("The Lorentz component is physical only if the Gaussian broadening model is justified.")
    if fwhm is None:
        warnings.append("Selected feature lacks two half-contrast crossings inside the fit window.")
    if selected == "peak" and q3 is None:
        warnings.append("No two absolute half-power crossings; Q_3dB is unavailable.")
    characteristic_width = min(width, fwhm) if fwhm else width
    local = np.abs(f - center) <= characteristic_width
    gaps = np.diff(f)
    # Count input-file intervals; an interpolated export does not establish
    # the solver's true frequency resolution.
    spacing = float(np.max(gaps[np.flatnonzero(local[:-1] | local[1:])])) if np.any(local) else float(np.max(gaps))
    samples_per_width = characteristic_width / spacing
    if samples_per_width < 8:
        warnings.append("Undersampled input: fewer than 8 input intervals across the relevant linewidth.")
    if min(center - f[0], f[-1] - center) < width:
        warnings.append("Resonance center is within one linewidth of a window edge; widen the fitting window.")
    if np.sqrt(rss / n) > .05 * max(abs(strength), np.finfo(float).eps):
        warnings.append("Residual RMSE exceeds 5% of fitted feature contrast; inspect multimode/background mismatch.")
    if np.any(fitted < -1e-6): warnings.append("Model predicts negative power inside the window.")
    hit_bound = bool(np.any(np.isclose(opt.x, lo, atol=1e-5)) or np.any(np.isclose(opt.x, hi, atol=1e-5)))
    if hit_bound: warnings.append("A shape parameter reached its fitting bound; the model may be unidentifiable.")
    diagnostics = {"n_points": n, "n_parameters": k, "rmse_power": float(np.sqrt(rss / n)),
                   "r_squared_power": 1 - rss / ss if ss > 0 else None,
                   "aicc": aicc, "bic": float(n * np.log(variance) + k * np.log(n)),
                   "optimizer_converged": bool(opt.success), "shape_at_bound": hit_bound,
                   "input_intervals_per_relevant_linewidth": float(samples_per_width),
                   "half_contrast_left": left, "half_contrast_right": right,
                   "feature_power": float(yd[idx]), "feature_contrast": strength,
                   "frequency_3db_peak": peak3,
                   "window": [float(f[0]), float(f[-1])]}
    return FitResult(model, frequency_unit, float(center), float(width), q_resonance,
                     None if q_spectral is None else float(q_spectral), ff, selected,
                     None if fwhm is None else float(fwhm),
                     None if q3 is None else float(q3), None if bw3 is None else float(bw3),
                     params, diagnostics, warnings, frequency=f, power=power,
                     fitted_power=fitted, fitted_baseline=base)


def fit_resonance(frequency, values, *, kind="amplitude", model="auto", window=None,
                  baseline="linear", feature="auto", frequency_unit="THz", n_starts=12):
    """Fit one isolated resonance. Inputs retain their frequency unit.

    auto ranks the five supported models by AICc, not by physical plausibility.
    Use an explicit window to isolate one mode from a multi-resonance spectrum.
    """
    if baseline not in ("constant", "linear") or feature not in ("auto", "peak", "dip"):
        raise ValueError("Invalid baseline or feature option.")
    if model != "auto" and model not in MODELS: raise ValueError(f"Unknown model: {model}")
    f = np.asarray(frequency, float)
    y = np.asarray(values)
    if f.ndim != 1 or y.ndim != 1 or len(f) != len(y):
        raise ValueError("frequency and values must be equal-length 1D arrays")
    if not np.all(np.isfinite(f)) or np.any(f <= 0):
        raise ValueError("Frequency must be finite and positive.")
    if not np.all(np.isfinite(y)): raise ValueError("Spectrum contains NaN or Inf.")
    if len(np.unique(f)) != len(f):
        raise ValueError("Duplicate frequencies: select one CST sweep block; do not average cases.")
    order = np.argsort(f); f, y = f[order], y[order]
    p = to_power(y, kind)
    if window is not None:
        if len(window) != 2 or window[0] >= window[1]: raise ValueError("window must be (low, high)")
        keep = (f >= window[0]) & (f <= window[1]); f, p = f[keep], p[keep]
    if len(f) < 12: raise ValueError("Need at least 12 samples in the fit window.")
    if np.ptp(p) <= max(1., float(np.max(np.abs(p)))) * 1e-12:
        raise ValueError("Spectrum is effectively flat; no linewidth can be estimated.")
    if baseline == "linear":
        xn = (f - np.mean(f)) / np.ptp(f)
        residual = p - np.polyval(np.polyfit(xn, p, 1), xn)
        if np.ptp(residual) <= max(1., float(np.max(np.abs(p)))) * 1e-12:
            raise ValueError("Spectrum is only a linear baseline; no resonance to fit.")
    models = MODELS if model == "auto" else (model,)
    results = [_one_fit(f, p, m, baseline, feature, frequency_unit, n_starts) for m in models]
    best = min(results, key=lambda r: r.diagnostics["aicc"] if r.diagnostics["aicc"] is not None else np.inf)
    if model == "auto":
        best.candidates = [{"model": r.model, "aicc": r.diagnostics["aicc"],
                            "rmse_power": r.diagnostics["rmse_power"],
                            "q_resonance": r.q_resonance, "q_spectral": r.q_spectral}
                           for r in sorted(results, key=lambda r: r.diagnostics["aicc"])]
        best.warnings.append("Automatic AICc ranking does not identify the physical resonance mechanism.")
    return best


def save_result(result, prefix, *, plot=True):
    """Save JSON, fitted sample table, fit/residual PNG. Never overwrite input."""
    prefix = Path(prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(str(prefix) + '.json').write_text(json.dumps(result.to_dict(), ensure_ascii=False,
                                                     indent=2, allow_nan=False), encoding='utf-8')
    with Path(str(prefix) + '.csv').open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.writer(handle)
        writer.writerow(['frequency_' + result.frequency_unit, 'power', 'fitted_power', 'baseline', 'residual'])
        writer.writerows(zip(result.frequency, result.power, result.fitted_power,
                             result.fitted_baseline, result.power - result.fitted_power))
    if not plot:
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    axes[0].plot(result.frequency, result.power, '.', ms=3, label='Input (power)')
    axes[0].plot(result.frequency, result.fitted_power, '-', label=result.model)
    axes[0].plot(result.frequency, result.fitted_baseline, '--', lw=1, label='Local baseline')
    axes[0].axvline(result.f0, color='gray', lw=.8, alpha=.7)
    axes[0].set_ylabel('Power response'); axes[0].legend()
    fmt = lambda value: 'N/A' if value is None else f'{value:.5g}'
    axes[0].set_title(f'{result.model}: pole Q={fmt(result.q_resonance)}, spectral Q={fmt(result.q_spectral)}')
    axes[1].plot(result.frequency, result.power - result.fitted_power, '.', ms=3)
    axes[1].axhline(0, color='gray', lw=.8)
    axes[1].set(xlabel=f'Frequency / {result.frequency_unit}', ylabel='Residual')
    fig.tight_layout(); fig.savefig(str(prefix) + '.png', dpi=180); plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='CSV, TXT or TSV frequency response')
    parser.add_argument('--kind', choices=['amplitude', 'power', 'db', 'complex'], default='amplitude')
    parser.add_argument('--model', choices=['auto', *MODELS], default='auto')
    parser.add_argument('--window', type=float, nargs=2, metavar=('LOW', 'HIGH'))
    parser.add_argument('--feature', choices=['auto', 'peak', 'dip'], default='auto')
    parser.add_argument('--baseline', choices=['constant', 'linear'], default='linear')
    parser.add_argument('--frequency-unit', default='THz')
    parser.add_argument('--fcol', type=int, default=0)
    parser.add_argument('--ycol', type=int, default=1)
    parser.add_argument('--imag-col', type=int)
    parser.add_argument('--skiprows', type=int, default=0)
    parser.add_argument('--delimiter')
    parser.add_argument('--starts', type=int, default=12)
    parser.add_argument('--output', default='results/qfit')
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args(argv)
    if (args.kind == 'complex') != (args.imag_col is not None):
        parser.error('--kind complex and --imag-col must be specified together')
    try:
        f, y = load_spectrum(args.input, args.fcol, args.ycol, args.imag_col, args.delimiter, args.skiprows)
        result = fit_resonance(f, y, kind=args.kind, model=args.model, window=args.window,
                               baseline=args.baseline, feature=args.feature,
                               frequency_unit=args.frequency_unit, n_starts=args.starts)
        source_path = Path(args.input).resolve()
        if any(Path(str(args.output) + ext).resolve() == source_path for ext in ('.csv', '.json', '.png')):
            raise ValueError('Output prefix would overwrite the input file.')
        save_result(result, args.output, plot=not args.no_plot)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(json.dumps(result.to_dict(), ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

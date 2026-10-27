# -*- coding: utf-8 -*-
from __future__ import annotations

"""Single-resonance power-spectrum fits with explicit linewidth/Q conventions.

Requires numpy, scipy and (only for plots) matplotlib. No signal is smoothed.
Fano uses eps=2*(f-f0)/Gamma, so the fitted pole Q is f0/Gamma.
"""

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

FIT_MODELS = ("lorentzian", "fano", "gaussian", "voigt", "pseudo_voigt")


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
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except UnicodeDecodeError:
        text = Path(path).read_text(encoding="gb18030")
    for lineno, line in enumerate(text.splitlines(), 1):
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
    if model != "auto" and model not in FIT_MODELS: raise ValueError(f"Unknown model: {model}")
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
    models = FIT_MODELS if model == "auto" else (model,)
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
    parser.add_argument('--model', choices=['auto', *FIT_MODELS], default='auto')
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


"""Multi-peak resonance measurements. Frequency: THz; input: |S21| amplitude.

half_prominence: power peak width at half its prominence (not half amplitude).
fano_pole: local power fit, epsilon=2*(f-f0)/Gamma, Q=f0/Gamma.
These two Q definitions need not agree for an asymmetric Fano resonance.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.optimize import least_squares


@dataclass
class Peak:
    index: int
    peak_frequency_THz: float
    peak_power: float
    prominence: float
    width_THz: float
    Q: float
    pole_frequency_THz: float = float("nan")
    fit_relative_rmse: float = float("nan")
    training_intervals_per_width: float = 0.0
    valid: bool = True
    reason: str = ""

    def to_dict(self):
        return asdict(self)


def _fit_fano(f, power, peak_index, left, right, initial_width):
    """Variable projection: fit center and positive linewidth, solve linear terms.

    P=b0+b1*x+C/(1+eps^2)+D*eps/(1+eps^2).
    Equivalent to a Fano pole plus a locally linear background; returns no q
    parameter because the amplitude/background decomposition is not unique.
    """
    xx, yy = f[left:right + 1], power[left:right + 1]
    if len(xx) < 12:
        return None
    step, span = float(np.max(np.diff(xx))), float(np.ptp(xx))
    origin = float(f[peak_index])
    xn = (xx - origin) / span
    gamma_lo, gamma_hi = max(step * 1.1, 1e-9), span * 1.5
    contrast = max(float(np.ptp(yy)), 1e-9)

    def project(z):
        center, gamma = z[0], np.exp(z[1])
        eps = 2 * (xx - center) / gamma
        mat = np.column_stack([np.ones_like(xx), xn, 1 / (1 + eps**2),
                               eps / (1 + eps**2)])
        coef = np.linalg.lstsq(mat, yy, rcond=None)[0]
        return mat @ coef

    candidates = []
    for factor, offset in [(1, 0), (0.6, -0.3), (1.6, 0.3)]:
        gamma = np.clip(initial_width * factor, gamma_lo * 1.01, gamma_hi * .99)
        center = np.clip(origin + offset * initial_width, xx[1], xx[-2])
        result = least_squares(lambda z: (project(z) - yy) / contrast,
                               [center, np.log(gamma)],
                               bounds=([xx[0], np.log(gamma_lo)],
                                       [xx[-1], np.log(gamma_hi)]),
                               max_nfev=70, ftol=1e-6, xtol=1e-7, gtol=1e-6)
        candidates.append((float(np.mean(result.fun**2)), result))
    err, result = min(candidates, key=lambda t: t[0])
    center, gamma = float(result.x[0]), float(np.exp(result.x[1]))
    # The fit must contain both sides of the pole, without width-bound pinning.
    ok = (result.success and center - xx[0] >= gamma * .6
          and xx[-1] - center >= gamma * .6
          and gamma > gamma_lo * 1.03 and gamma < gamma_hi * .97)
    return center, gamma, float(np.sqrt(err)), ok


def extract_peaks(frequency, amplitude, q_definition="fano_pole", prominence=.01,
                  min_power=.05, training_step=None, min_intervals=8,
                  fit_max_error=.12):
    f, a = np.asarray(frequency, float), np.asarray(amplitude, float)
    if (f.ndim != 1 or a.shape != f.shape or len(f) < 5
            or not np.isfinite(f).all() or not np.isfinite(a).all()
            or np.any(np.diff(f) <= 0) or np.any(a < 0)):
        raise ValueError("频率必须递增且谱线有限、非负；需要至少五个频点。")
    if q_definition not in ("fano_pole", "half_prominence"):
        raise ValueError("未知 Q 定义")
    power = a**2
    indices, props = find_peaks(power, prominence=prominence, height=min_power)
    if not len(indices):
        return []
    widths = peak_widths(power, indices, rel_height=.5,
                        prominence_data=(props["prominences"],
                                         props["left_bases"], props["right_bases"]))
    step = float(training_step if training_step is not None else np.max(np.diff(f)))
    if step <= 0:
        raise ValueError("训练频率步长必须大于 0。")
    measured = []
    for j, idx in enumerate(indices):
        lpos, rpos = float(widths[2][j]), float(widths[3][j])
        lf = float(np.interp(lpos, np.arange(len(f)), f))
        rf = float(np.interp(rpos, np.arange(len(f)), f))
        width = rf - lf
        peak = Peak(int(idx), float(f[idx]), float(power[idx]),
                    float(props["prominences"][j]), width, float(f[idx] / width))
        reasons = []
        # A truncated peak/insufficient flanks must not produce a confident Q.
        if lpos < 1 or rpos > len(f) - 2:
            reasons.append("peak_at_window_edge")
        if q_definition == "fano_pole":
            span = max(4 * width, 12 * float(np.median(np.diff(f))))
            left = int(np.searchsorted(f, f[idx] - span))
            right = min(len(f) - 1, int(np.searchsorted(f, f[idx] + span)))
            if j:
                left = max(left, int((indices[j-1] + idx) // 2))
            if j + 1 < len(indices):
                right = min(right, int((idx + indices[j+1]) // 2))
            fit = _fit_fano(f, power, int(idx), left, right, width)
            if fit is None:
                reasons.append("too_few_fit_points")
            else:
                f0, gamma, err, fit_ok = fit
                peak.pole_frequency_THz, peak.width_THz = f0, gamma
                peak.fit_relative_rmse, peak.Q = err, f0 / gamma
                if not fit_ok or err > fit_max_error:
                    reasons.append("local_fano_fit_unreliable_or_overlapping")
        peak.training_intervals_per_width = peak.width_THz / step
        if peak.training_intervals_per_width < min_intervals:
            reasons.append("training_grid_under_resolves_width")
        peak.valid = not reasons
        peak.reason = ";".join(reasons)
        measured.append(peak)
    return measured


def match_peaks(peaks, target_frequency, target_q, frequency_tolerance,
                q_tolerance, policy="penalize", extra_relative=.5,
                extra_weight=.2, frequency_reference="peak"):
    """Compare every peak against the target, not just the largest peak.

    Returns score, selected peak, extra count, and whether hard criteria pass.
    Invalid peaks retain a finite penalty to guide the search, but cannot pass.
    """
    if not peaks:
        return 1e6, None, 0, False
    if policy not in ("match", "penalize", "single", "dominant"):
        raise ValueError("未知多峰策略")
    candidates = []
    for peak in peaks:
        freq = (peak.peak_frequency_THz if frequency_reference == "peak"
                else peak.pole_frequency_THz)
        if not np.isfinite(freq) or not np.isfinite(peak.Q):
            continue
        ef = abs(freq - target_frequency) / frequency_tolerance
        eq = abs(peak.Q - target_q) / q_tolerance
        extras = sum(p.index != peak.index and p.prominence >= peak.prominence * extra_relative
                     for p in peaks)
        dominant = peak.peak_power >= max(p.peak_power for p in peaks) - 1e-10
        hard_ok = peak.valid and (policy != "single" or extras == 0)
        hard_ok = hard_ok and (policy != "dominant" or dominant)
        score = float(ef**2 + eq**2)
        if policy == "penalize":
            score += extra_weight * extras
        if not hard_ok:
            score += 100 + (100 * extras if policy == "single" else 0)
        candidates.append((score, peak, int(extras), bool(hard_ok and ef <= 1 and eq <= 1)))
    return min(candidates, key=lambda t: t[0]) if candidates else (1e6, None, 0, False)


"""Power-spectrum Q measurement, reusing the forward-screening definitions."""
from dataclasses import dataclass, field, asdict
import numpy as np
from scipy.signal import find_peaks, peak_widths

MODELS = {"法诺 Fano": "fano", "洛伦兹 Lorentzian": "lorentzian",
          "高斯 Gaussian": "gaussian", "Voigt": "voigt",
          "伪 Voigt": "pseudo_voigt", "自动比较 AICc": "auto"}
METHODS = ["法诺极点 Q（筛选同算法）", "半突出度 Q（筛选同算法）", "自定义线形拟合"]
KINDS = {"线性幅度 |S|": "amplitude", "功率 |S|²": "power",
         "分贝 dB": "db", "复数：实部+虚部": "complex"}
UNITS = {"THz": 1., "GHz": .001, "MHz": 1e-6, "Hz": 1e-12}


def read_data(path, kind="amplitude", unit="THz", fcol=1, ycol=2, imagcol=3):
    if min(fcol,ycol,imagcol) < 1:
        raise ValueError("列号从 1 开始。")
    f,y = load_spectrum(path, fcol=fcol-1,ycol=ycol-1,
                        imag_col=imagcol-1 if kind=="complex" else None)
    f=f*UNITS[unit]
    if not np.isfinite(f).all() or np.any(f<=0) or len(f)<5:
        raise ValueError("频率需为正数，且至少有五个有效频点。")
    if len(np.unique(f)) != len(f):
        raise ValueError("频率重复：文件可能包含多个扫参数据块，请先导出单组曲线。")
    order=np.argsort(f)
    return f[order], to_power(y[order],kind)


def detect_features(f,power,prominence=.01,min_power=.05,feature="peak"):
    if prominence <= 0 or min_power < 0:
        raise ValueError("突出度须大于零，最低功率不能为负数。")
    result=[]
    for direction in (["peak","dip"] if feature=="both" else [feature]):
        sign=1 if direction=="peak" else -1
        indices,props=find_peaks(sign*power,prominence=prominence)
        if not len(indices): continue
        widths=peak_widths(sign*power,indices,rel_height=.5,
            prominence_data=(props['prominences'],props['left_bases'],props['right_bases']))
        for j,idx in enumerate(indices):
            if direction=='peak' and power[idx]<min_power: continue
            left=float(np.interp(widths[2][j],np.arange(len(f)),f))
            right=float(np.interp(widths[3][j],np.arange(len(f)),f))
            result.append(dict(index=int(idx),frequency=float(f[idx]),power=float(power[idx]),
                prominence=float(props['prominences'][j]),kind=direction,left=left,right=right,
                half_level=float(sign*widths[1][j]),width=right-left))
    return sorted(result,key=lambda v:v['frequency'])


def local_window(f,power,feature,prominence=.01,min_power=.05):
    """Same local window rule used by extract_peaks for a transmission peak."""
    all_features=detect_features(f,power,prominence,min_power,feature['kind'])
    j=min(range(len(all_features)),key=lambda k:abs(all_features[k]['frequency']-feature['frequency']))
    idx=all_features[j]['index']
    width=all_features[j]['width']
    span=max(4*width,12*float(np.median(np.diff(f))))
    left=int(np.searchsorted(f,f[idx]-span))
    right=min(len(f)-1,int(np.searchsorted(f,f[idx]+span)))
    if j: left=max(left,(all_features[j-1]['index']+idx)//2)
    if j+1<len(all_features): right=min(right,(idx+all_features[j+1]['index'])//2)
    return left,right


@dataclass
class Measurement:
    method: str
    model: str
    primary_q: float | None
    q_definition: str
    center: float
    visible_frequency: float
    linewidth: float
    q_pole: float | None = None
    q_half_contrast: float | None = None
    q_3db: float | None = None
    q_half_prominence: float | None = None
    diagnostics: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    f: np.ndarray = field(default_factory=lambda: np.array([]))
    power: np.ndarray = field(default_factory=lambda: np.array([]))
    fitted: np.ndarray = field(default_factory=lambda: np.array([]))
    baseline: np.ndarray = field(default_factory=lambda: np.array([]))
    guide: dict = field(default_factory=dict)
    details: dict = field(default_factory=dict)

    def to_dict(self):
        result=asdict(self)
        for key in ('f','power','fitted','baseline'): result.pop(key)
        return result


def measure(f,power,method=METHODS[0],model='fano',baseline='linear',feature='peak',
            target_frequency=None,prominence=.01,min_power=.05):
    f,power=np.asarray(f,float),np.asarray(power,float)
    if (f.ndim!=1 or power.shape!=f.shape or len(f)<5 or np.any(np.diff(f)<=0)
        or not np.isfinite(f).all() or not np.isfinite(power).all() or np.any(power<0)):
        raise ValueError("频率需递增，功率数据需非负且有限，窗口内至少五个点。")
    features=detect_features(f,power,prominence,min_power,feature)
    target=(min(features,key=lambda v:abs(v['frequency']-target_frequency))
            if target_frequency is not None and features else
            max(features,key=lambda v:v['prominence']) if features else None)
    if method in METHODS[:2]:
        if not target: raise ValueError("窗口内未检测到目标峰/谷。请检查窗口或降低突出度/最低功率。")
        if method==METHODS[0]:
            if feature!='peak':
                raise ValueError("筛选脚本的法诺同算法用于透射峰。拟合透射谷请选择‘自定义线形拟合’。")
            peaks=extract_peaks(f,np.sqrt(power),'fano_pole',prominence,min_power,
                                training_step=float(np.max(np.diff(f))))
            peak=min(peaks,key=lambda pk:abs(pk.peak_frequency_THz-target['frequency']))
            if not np.isfinite(peak.pole_frequency_THz):
                raise ValueError("窗口内没有足够点进行法诺拟合。请扩大窗口或加密真实采样。")
            left,right=local_window(f,power,target,prominence,min_power)
            xx,yy=f[left:right+1],power[left:right+1]
            xn=(xx-target['frequency'])/np.ptp(xx)
            eps=2*(xx-peak.pole_frequency_THz)/peak.width_THz
            design=np.column_stack([np.ones_like(xx),xn,1/(1+eps**2),eps/(1+eps**2)])
            coefficients=np.linalg.lstsq(design,yy,rcond=None)[0]
            fitted=design@coefficients
            messages=[]
            relative_feature_rmse=float(np.sqrt(np.mean((fitted-yy)**2))/max(peak.prominence,1e-12))
            if relative_feature_rmse>.05:
                messages.append('拟合误差相对峰突出度较大，请调整窗口并检查 Q 是否稳定；背景可能非线性。')
            reasons={'training_grid_under_resolves_width':'一个线宽内不足 8 个输入频率间隔，Q 的可靠性不足。',
                'local_fano_fit_unreliable_or_overlapping':'局部法诺拟合不可靠，可能有峰重叠、边界不足或线形不符。',
                'peak_at_window_edge':'峰靠近窗口边缘，请扩大窗口。',
                'too_few_fit_points':'拟合频点不足。'}
            for code in peak.reason.split(';'):
                if code: messages.append(reasons.get(code,code))
            return Measurement(method,'fano',peak.Q,'极点 Q=f0/Gamma',peak.pole_frequency_THz,
                peak.peak_frequency_THz,peak.width_THz,q_pole=peak.Q,
                diagnostics=dict(valid=peak.valid,relative_rmse=peak.fit_relative_rmse,
                    rmse_relative_to_peak_prominence=relative_feature_rmse,
                    intervals_per_width=peak.training_intervals_per_width,
                    fit_window=[float(xx[0]),float(xx[-1])]),warnings=messages,
                f=xx,power=yy,fitted=fitted,baseline=design[:,:2]@coefficients[:2],
                details=dict(peak=peak.to_dict(),coefficients=coefficients.tolist(),
                    epsilon='2*(f-f0)/Gamma',x_origin=target['frequency'],x_scale=float(np.ptp(xx))))
        width=target['width']; q=target['frequency']/width
        intervals=width/np.max(np.diff(f)); messages=[]
        if intervals<8: messages.append('一个线宽内不足 8 个输入频率间隔，Q 的可靠性不足。')
        if target['left']<=f[1] or target['right']>=f[-2]:
            messages.append('半突出度交点靠近边界，请扩大窗口复核。')
        if target['kind']=='dip':
            messages.append('当前测量透射谷的半突出度宽度；原筛选脚本只检测透射峰。')
        return Measurement(method,'no_fit',q,'半突出度 Q=f_feature/width',target['frequency'],
            target['frequency'],width,q_half_prominence=q,
            diagnostics=dict(intervals_per_width=float(intervals)),warnings=messages,
            f=f,power=power,guide=dict(left=target['left'],right=target['right'],level=target['half_level']))
    if method!=METHODS[2]: raise ValueError('未知计算方式。')
    result=fit_resonance(f,power,kind='power',model=model,baseline=baseline,feature=feature,
                         frequency_unit='THz',n_starts=8)
    warnings=[]
    if len(features)>1:
        warnings.append('当前区间有多个显著峰/谷，单共振拟合可能无效；建议用“选中峰邻域”缩小区间。')
    if target_frequency is not None and target and abs(result.feature_frequency-target_frequency)>max(target['width'],2*np.max(np.diff(f))):
        warnings.append('拟合特征偏离你选中的峰/谷，请缩小拟合窗口。')
    warning_map={
        'Optimizer did not converge.':'优化未收敛。',
        'Spectral Q is a linewidth descriptor; it is not automatically a cavity decay Q.':'高斯/Voigt 等得到的是谱宽 Q，不自动等于模态衰减 Q。',
        'Fano pole Q=f0/Gamma and feature-width Q have different definitions.':'Fano 的极点 Q 与峰形宽度 Q 定义不同。',
        'The Lorentz component is physical only if the Gaussian broadening model is justified.':'Voigt 中洛伦兹分量的物理解释需有展宽机制依据。',
        'Selected feature lacks two half-contrast crossings inside the fit window.':'窗口内缺少完整的半对比度交点，该 Q 不可用。',
        'No two absolute half-power crossings; Q_3dB is unavailable.':'窗口内没有两个绝对半功率交点，3dB Q 不可用。',
        'Undersampled input: fewer than 8 input intervals across the relevant linewidth.':'相关线宽内不足 8 个输入频率间隔，采样分辨率不足。',
        'Resonance center is within one linewidth of a window edge; widen the fitting window.':'共振中心距离区间边缘不足一个线宽，建议扩大窗口。',
        'Residual RMSE exceeds 5% of fitted feature contrast; inspect multimode/background mismatch.':'拟合误差大于特征对比度的 5%，请检查多峰或背景不符。',
        'Model predicts negative power inside the window.':'拟合出现负功率，线形/窗口可能不合适。',
        'A shape parameter reached its fitting bound; the model may be unidentifiable.':'参数碰到拟合边界，Q 可能不稳定。',
        'Automatic AICc ranking does not identify the physical resonance mechanism.':'自动 AICc 只比较拟合，不判断实际共振机理。'}
    warnings.extend(warning_map.get(w,w) for w in result.warnings)
    primary=result.q_resonance if result.q_resonance is not None else result.q_spectral
    return Measurement(method,result.model,primary,
        '极点 Q=f0/Gamma' if result.q_resonance is not None else '拟合半对比度谱宽 Q',
        result.f0,result.feature_frequency,result.linewidth,q_pole=result.q_resonance,
        q_half_contrast=result.q_spectral,q_3db=result.q_3db,diagnostics=result.diagnostics,
        warnings=warnings,f=result.frequency,power=result.power,fitted=result.fitted_power,
        baseline=result.fitted_baseline,details=result.to_dict(),
        guide=dict(left=result.diagnostics['half_contrast_left'],
                   right=result.diagnostics['half_contrast_right']))


def clean_json(value):
    if isinstance(value,dict): return {k:clean_json(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,np.ndarray)): return [clean_json(v) for v in value]
    if isinstance(value,(float,np.floating)): return float(value) if np.isfinite(value) else None
    if isinstance(value,np.integer): return int(value)
    if isinstance(value,np.bool_): return bool(value)
    return value


# -*- coding: utf-8 -*-
"""Interactive resonance Q calculator. Run: python qcalc_gui.py [spectrum.txt]."""
import argparse
import csv
import json
import queue
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import SpanSelector

ROOT=Path(__file__).resolve().parent
matplotlib.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus']=False


def fmt(value):
    return '不可用' if value is None or not np.isfinite(value) else f'{value:.7g}'


class QCalculator:
    def __init__(self,root):
        self.root=root; self.root.title('谐振 Q 值计算器 · 多峰选取与线形拟合')
        self.root.geometry('1420x940'); self.root.minsize(1120,760)
        self.path=None; self.f=None; self.power=None; self.features=[]
        self.result=None; self.snapshot=None; self.revision=0; self.busy=False
        self.jobs=queue.Queue(); self.loaded_options=None; self.highlight=None
        self.selection_signature=None; self.selector=None
        self.filename=tk.StringVar(value='打开 TXT / CSV 数据，或先加载示例')
        self.kind=tk.StringVar(value=next(iter(KINDS)))
        self.unit=tk.StringVar(value='THz')
        self.fcol=tk.StringVar(value='1'); self.ycol=tk.StringVar(value='2'); self.icol=tk.StringVar(value='3')
        self.low=tk.StringVar(value='0.1'); self.high=tk.StringVar(value='0.5')
        self.method=tk.StringVar(value=METHODS[0]); self.model=tk.StringVar(value='法诺 Fano')
        self.baseline=tk.StringVar(value='线性背景'); self.feature_kind=tk.StringVar(value='透射峰')
        self.prominence=tk.StringVar(value='0.01'); self.min_power=tk.StringVar(value='0.05')
        self.status=tk.StringVar(value='就绪：先打开数据。内部统一用功率谱拟合，频率窗口单位为 THz。')
        self.q_text=tk.StringVar(value='Q = —')
        self._build()
        for variable in (self.kind,self.unit,self.fcol,self.ycol,self.icol,self.low,self.high,
                         self.method,self.model,self.baseline,self.feature_kind,self.prominence,self.min_power):
            variable.trace_add('write',self.invalidate)
        self.method.trace_add('write',self.update_modes)
        self.update_modes(); self.root.after(80,self.poll)

    def _build(self):
        style=ttk.Style(); style.theme_use('clam')
        style.configure('TButton',padding=(8,5)); style.configure('TLabel',font=('Microsoft YaHei UI',10))
        style.configure('TLabelframe.Label',font=('Microsoft YaHei UI',10,'bold'))
        top=ttk.Frame(self.root,padding=(12,8)); top.pack(fill='x')
        ttk.Button(top,text='打开数据…',command=self.open_dialog).pack(side='left')
        ttk.Button(top,text='加载演示谱',command=self.load_demo).pack(side='left',padx=8)
        ttk.Label(top,textvariable=self.filename).pack(side='left',fill='x',expand=True,padx=5)
        ttk.Label(self.root,textvariable=self.status,wraplength=1350,padding=(12,8),foreground='#37485f').pack(side='bottom',fill='x')
        main=ttk.Panedwindow(self.root,orient='horizontal'); main.pack(fill='both',expand=True,padx=12)
        left=ttk.Frame(main,width=350); right=ttk.Frame(main)
        main.add(left,weight=0); main.add(right,weight=1)
        self.root.after(40,lambda:main.sashpos(0,420))
        data=ttk.LabelFrame(left,text='1  数据设置',padding=9); data.pack(fill='x',pady=(0,8))
        self.combo(data,'数据类型',self.kind,list(KINDS),0)
        self.combo(data,'原文件频率单位',self.unit,list(UNITS),1)
        for row,(label,var) in enumerate([('频率列（从 1 开始）',self.fcol),('数值 / 实部列',self.ycol),('虚部列（仅复数）',self.icol)],2):
            ttk.Label(data,text=label).grid(row=row,column=0,sticky='w',pady=3)
            ttk.Entry(data,textvariable=var,width=16).grid(row=row,column=1,sticky='ew',padx=(8,0))
        ttk.Button(data,text='应用格式 / 重新读取',command=self.reload).grid(row=5,column=0,columnspan=2,sticky='ew',pady=(5,0))
        win=ttk.LabelFrame(left,text='2  频率窗口（THz）',padding=9); win.pack(fill='x',pady=(0,8))
        for row,(label,var) in enumerate([('起始频率',self.low),('结束频率',self.high)]):
            ttk.Label(win,text=label).grid(row=row,column=0,sticky='w',pady=3)
            ttk.Entry(win,textvariable=var,width=18).grid(row=row,column=1,sticky='ew',padx=8)
        btns=ttk.Frame(win); btns.grid(row=2,column=0,columnspan=2,sticky='ew',pady=4)
        ttk.Button(btns,text='应用窗口',command=self.apply_window).pack(side='left')
        ttk.Button(btns,text='全频段',command=self.full_window).pack(side='left',padx=4)
        ttk.Label(win,text='也可在主图按住左键拖动选区。',foreground='#566477').grid(row=3,column=0,columnspan=2,sticky='w')
        calc=ttk.LabelFrame(left,text='3  计算方式与拟合线形',padding=9); calc.pack(fill='x',pady=(0,8))
        self.combo(calc,'计算方式',self.method,METHODS,0,width=26)
        self.model_box=self.combo(calc,'拟合线形',self.model,list(MODELS),1)
        self.base_box=self.combo(calc,'拟合背景',self.baseline,['线性背景','常数背景'],2)
        self.combo(calc,'目标特征',self.feature_kind,['透射峰','透射谷'],3)
        for row,(label,var) in enumerate([('最小突出度（功率）',self.prominence),('最低峰功率',self.min_power)],4):
            ttk.Label(calc,text=label).grid(row=row,column=0,sticky='w',pady=3)
            ttk.Entry(calc,textvariable=var,width=14).grid(row=row,column=1,sticky='ew',padx=(8,0))
        actions=ttk.Frame(calc); actions.grid(row=6,column=0,columnspan=2,sticky='ew',pady=(6,0))
        self.fit_button=ttk.Button(actions,text='计算 Q / 拟合',command=self.start_fit)
        self.fit_button.pack(side='left',fill='x',expand=True)
        self.export_button=ttk.Button(actions,text='导出结果…',command=self.export,state='disabled')
        self.export_button.pack(side='left',padx=(6,0))
        res=ttk.LabelFrame(left,text='计算结果',padding=9); res.pack(fill='both',expand=True)
        ttk.Label(res,textvariable=self.q_text,font=('Microsoft YaHei UI',20,'bold'),foreground='#145caa').pack(anchor='w')
        result_frame=ttk.Frame(res); result_frame.pack(fill='both',expand=True,pady=(6,0))
        self.result_box=tk.Text(result_frame,height=9,width=30,wrap='word',font=('Microsoft YaHei UI',10),relief='flat',bg='#f4f7fb')
        result_scroll=ttk.Scrollbar(result_frame,orient='vertical',command=self.result_box.yview)
        self.result_box.configure(yscrollcommand=result_scroll.set)
        result_scroll.pack(side='right',fill='y'); self.result_box.pack(fill='both',expand=True)
        self.write_result('主结果随计算方式而定。\nFano 极点 Q 与谱形宽度 Q 不等价。')
        self.fig=Figure(figsize=(9,5.8),dpi=100,layout='constrained')
        grid=self.fig.add_gridspec(2,1,height_ratios=[3,1])
        self.ax=self.fig.add_subplot(grid[0]); self.resax=self.fig.add_subplot(grid[1],sharex=self.ax)
        self.canvas=FigureCanvasTkAgg(self.fig,master=right); self.canvas.get_tk_widget().pack(fill='both',expand=True)
        toolbar=NavigationToolbar2Tk(self.canvas,right,pack_toolbar=False); toolbar.update(); toolbar.pack(fill='x')
        self.toolbar=toolbar
        bottom=ttk.LabelFrame(right,text='窗口内的峰 / 谷：点选目标后计算；可缩到该峰邻域',padding=7)
        bottom.pack(fill='x',pady=(5,0))
        header=ttk.Frame(bottom); header.pack(fill='x')
        ttk.Button(header,text='刷新检测',command=self.refresh_features).pack(side='left')
        ttk.Button(header,text='选中峰邻域',command=self.isolate_selected).pack(side='left',padx=6)
        self.feature_count=ttk.Label(header,text='尚未加载数据'); self.feature_count.pack(side='left',padx=6)
        table_frame=ttk.Frame(bottom); table_frame.pack(fill='x',pady=(5,0))
        columns=('type','frequency','power','prominence','q_width')
        self.table=ttk.Treeview(table_frame,columns=columns,show='headings',height=4,selectmode='browse')
        for col,title in zip(columns,['类型','频率 / THz','功率','突出度','半突出度 Q']):
            self.table.heading(col,text=title); self.table.column(col,width=95,anchor='center')
        scroll=ttk.Scrollbar(table_frame,orient='vertical',command=self.table.yview)
        self.table.configure(yscrollcommand=scroll.set); scroll.pack(side='right',fill='y'); self.table.pack(fill='x',expand=True)
        self.table.bind('<<TreeviewSelect>>',self.feature_selected)
        self.draw()

    def combo(self,parent,label,variable,values,row,width=20):
        ttk.Label(parent,text=label).grid(row=row,column=0,sticky='w',pady=4)
        box=ttk.Combobox(parent,textvariable=variable,values=values,state='readonly',width=width)
        box.grid(row=row,column=1,sticky='ew',padx=(8,0)); parent.columnconfigure(1,weight=1)
        return box

    def invalidate(self,*_):
        self.revision+=1; self.result=None; self.snapshot=None
        self.export_button.configure(state='disabled'); self.q_text.set('Q = —（待计算）')
        self.write_result('设置已更改，请点击“计算 Q / 拟合”。')

    def update_modes(self,*_):
        custom=self.method.get()==METHODS[2]
        self.model_box.configure(state='readonly' if custom else 'disabled')
        self.base_box.configure(state='readonly' if custom else 'disabled')

    def data_options(self):
        return dict(kind=KINDS[self.kind.get()],unit=self.unit.get(),fcol=int(self.fcol.get()),
                    ycol=int(self.ycol.get()),imagcol=int(self.icol.get()))

    def feature_options(self):
        return dict(prominence=float(self.prominence.get()),min_power=float(self.min_power.get()),
                    feature='peak' if self.feature_kind.get()=='透射峰' else 'dip')

    def open_dialog(self):
        path=filedialog.askopenfilename(title='打开单组频谱数据',filetypes=[('频谱数据','*.txt *.csv *.tsv *.dat'),('所有文件','*.*')])
        if path: self.load_file(path)

    def load_file(self,path):
        try:
            options=self.data_options(); f,power=read_data(path,**options)
            self.path=Path(path); self.f,self.power=f,power; self.loaded_options=options
            self.filename.set(str(self.path)); self.full_window()
            self.status.set(f'已加载 {len(f)} 个点；数值按 {self.kind.get()} 转为功率；频率已统一为 THz。')
        except Exception as exc: self.error(exc)

    def reload(self):
        if self.path: self.load_file(self.path)
        else: self.open_dialog()

    def load_demo(self):
        f=np.linspace(.1,.5,2001); eps=2*(f-.30)/.008
        self.f=f; self.power=.05+.12*(1.8+eps)**2/(1+eps**2)+.22/(1+(2*(f-.42)/.016)**2)
        self.path=None; self.loaded_options=self.data_options()
        self.filename.set('内置演示：0.30 THz 附近法诺峰 + 0.42 THz 洛伦兹峰')
        self.full_window(); self.status.set('演示谱已加载。点选峰，或拖动图中区间，再点击计算。')

    def window_data(self):
        if self.f is None: raise ValueError('请先打开数据或加载示例。')
        lo,hi=float(self.low.get()),float(self.high.get())
        if not np.isfinite([lo,hi]).all() or lo>=hi: raise ValueError('起始频率必须小于结束频率。')
        mask=(self.f>=lo)&(self.f<=hi)
        if mask.sum()<5: raise ValueError('窗口内不足五个数据点，请扩大频率范围。')
        return self.f[mask],self.power[mask]

    def full_window(self):
        if self.f is None: return
        self.low.set(f'{self.f[0]:.10g}'); self.high.set(f'{self.f[-1]:.10g}')
        self.apply_window()

    def span_selected(self,low,high):
        if self.toolbar.mode or self.f is None or high-low<=0: return
        self.low.set(f'{low:.10g}'); self.high.set(f'{high:.10g}'); self.apply_window()

    def apply_window(self):
        try:
            self.window_data(); self.invalidate(); self.refresh_features()
        except Exception as exc: self.error(exc)

    def refresh_features(self):
        try:
            f,p=self.window_data(); opts=self.feature_options()
            features=detect_features(f,p,**opts)
            self.features=features
            self.table.delete(*self.table.get_children())
            for j,item in enumerate(features):
                self.table.insert('', 'end',iid=str(j),values=(
                    '峰' if item['kind']=='peak' else '谷',f"{item['frequency']:.7g}",
                    f"{item['power']:.5g}",f"{item['prominence']:.5g}",f"{item['frequency']/item['width']:.5g}"))
            if features:
                j=max(range(len(features)),key=lambda k:features[k]['prominence'])
                self.table.selection_set(str(j)); self.table.see(str(j))
            self.feature_count.configure(text=f'检测到 {len(features)} 个特征')
            self.draw()
        except Exception as exc: self.error(exc)

    def selected_feature(self):
        selection=self.table.selection()
        return self.features[int(selection[0])] if selection else None

    def feature_selected(self,event=None):
        selected=self.selected_feature()
        signature=(selected['kind'],selected['frequency']) if selected else None
        if signature!=self.selection_signature:
            self.selection_signature=signature
            self.invalidate(); self.draw()

    def isolate_selected(self):
        try:
            f,p=self.window_data(); item=self.selected_feature()
            if item is None: raise ValueError('请先在表格中选择一个峰或谷。')
            opts=self.feature_options()
            left,right=local_window(f,p,item,opts['prominence'],opts['min_power'])
            self.low.set(f'{f[left]:.10g}'); self.high.set(f'{f[right]:.10g}'); self.apply_window()
        except Exception as exc: self.error(exc)

    def start_fit(self):
        if self.busy: return
        try:
            if self.path and self.data_options()!=self.loaded_options:
                raise ValueError('数据格式已更改，请先点击“应用格式 / 重新读取”。')
            f,p=self.window_data(); opts=self.feature_options(); item=self.selected_feature()
            if not np.isfinite([opts['prominence'],opts['min_power']]).all():
                raise ValueError('突出度和最低峰功率必须是有限数值。')
            kwargs=dict(method=self.method.get(),model=MODELS[self.model.get()],
                baseline='linear' if self.baseline.get()=='线性背景' else 'constant',
                target_frequency=item['frequency'] if item else None,**opts)
            self.invalidate(); rev=self.revision
            snapshot=dict(source=str(self.path) if self.path else 'built-in demo',
                data_options=self.loaded_options,user_window_THz=[float(self.low.get()),float(self.high.get())],
                selected_feature=item,settings=kwargs)
            self.busy=True; self.fit_button.configure(state='disabled')
            self.status.set('正在计算… 窗口仍可操作；更改设置后旧结果会自动作废。')
            def worker():
                try: self.jobs.put((rev,measure(f.copy(),p.copy(),**kwargs),snapshot,None))
                except Exception as exc: self.jobs.put((rev,None,None,exc))
            threading.Thread(target=worker,daemon=True).start()
        except Exception as exc: self.error(exc)

    def poll(self):
        try:
            while True:
                rev,result,snapshot,error=self.jobs.get_nowait()
                self.busy=False; self.fit_button.configure(state='normal')
                if rev!=self.revision:
                    self.status.set('计算期间设置已改变，旧结果已丢弃；请按当前设置重新计算。'); continue
                if error: self.error(error); continue
                self.result=result; self.snapshot=snapshot; self.show_result()
        except queue.Empty: pass
        self.root.after(80,self.poll)

    def write_result(self,text):
        self.result_box.configure(state='normal'); self.result_box.delete('1.0','end')
        self.result_box.insert('1.0',text); self.result_box.configure(state='disabled')

    def show_result(self):
        r=self.result
        self.q_text.set('Q = '+fmt(r.primary_q))
        text=(f'定义：{r.q_definition}\n线形：{r.model}\n'
              f'中心 / 极点：{fmt(r.center)} THz\n可见特征：{fmt(r.visible_frequency)} THz\n'
              f'模型线宽：{fmt(r.linewidth*1000)} GHz\n\n'
              f'极点 Q：{fmt(r.q_pole)}\n半突出度 Q：{fmt(r.q_half_prominence)}\n'
              f'拟合半对比度 Q：{fmt(r.q_half_contrast)}\n绝对半功率 Q：{fmt(r.q_3db)}\n')
        if 'r_squared_power' in r.diagnostics: text+=f"功率拟合 R²：{fmt(r.diagnostics['r_squared_power'])}\n"
        if 'relative_rmse' in r.diagnostics: text+=f"相对拟合误差：{fmt(r.diagnostics['relative_rmse'])}\n"
        if r.warnings: text+='\n注意：\n'+'\n'.join('• '+w for w in r.warnings)
        self.write_result(text); self.export_button.configure(state='normal')
        self.status.set('计算完成。'+ ('存在可靠性提示，请查看结果栏。' if r.warnings else '可导出曲线、结果和图像。'))
        self.draw()

    def draw(self):
        self.ax.clear(); self.resax.clear()
        self.ax.set(ylabel='功率响应 |S|²',title='原始数据与拟合曲线')
        self.resax.set(xlabel='频率 / THz',ylabel='残差')
        if self.f is not None:
            self.ax.plot(self.f,self.power,color='#98a9bd',lw=1,label='原始功率谱')
            try:
                f,p=self.window_data(); self.ax.plot(f,p,'.',color='#175bbb',ms=2,label='选中窗口')
                lo,hi=float(self.low.get()),float(self.high.get())
                self.ax.set_xlim(lo,hi)
                pmin,pmax=float(p.min()),float(p.max()); margin=max((pmax-pmin)*.13,.005)
                self.ax.set_ylim(pmin-margin,pmax+margin)
            except ValueError: pass
            item=self.selected_feature()
            if item:
                self.ax.plot(item['frequency'],item['power'],'o',color='#e18223',ms=7,label='选中特征')
        if self.result is not None:
            r=self.result
            if len(r.fitted):
                self.ax.plot(r.f,r.fitted,color='#e04a46',lw=2,label=f'拟合：{r.model}')
                self.ax.plot(r.f,r.baseline,'--',color='#45a17f',lw=1,label='拟合背景')
                self.resax.plot(r.f,r.power-r.fitted,color='#6c5ca7',lw=1)
            if r.guide.get('left') is not None and r.guide.get('right') is not None:
                for x in (r.guide['left'],r.guide['right']): self.ax.axvline(x,color='#9b71b8',ls=':',lw=1)
                if 'level' in r.guide:
                    self.ax.hlines(r.guide['level'],r.guide['left'],r.guide['right'],color='#9b71b8',lw=2,label='半突出度宽度')
            self.ax.axvline(r.center,color='#616f80',ls=':',lw=1,label='中心 / 极点')
            self.ax.set_title(f'{r.model} · Q={fmt(r.primary_q)} · {r.q_definition}',fontsize=11)
        self.resax.axhline(0,color='#a7afbc',lw=.8)
        for ax in (self.ax,self.resax): ax.grid(alpha=.18)
        if self.f is not None: self.ax.legend(loc='best',fontsize=8)
        if self.selector is not None: self.selector.disconnect_events()
        self.selector=SpanSelector(self.ax,self.span_selected,'horizontal',useblit=True,
                                  props=dict(alpha=.17,facecolor='#2391d0'),interactive=True)
        self.canvas.draw_idle()

    def export(self,prefix=None,overwrite=False):
        if self.result is None:
            messagebox.showinfo('尚无结果','请先计算。'); return
        if prefix is None:
            stem=self.path.stem if self.path else 'demo'
            prefix=filedialog.asksaveasfilename(title='选择保存名称（自动生成 JSON、CSV 和 PNG）',
                initialfile=stem+'_q_result',defaultextension='.json',filetypes=[('结果文件','*.json')])
            if not prefix: return
        try:
            prefix=Path(prefix)
            if prefix.suffix.lower() in ('.json','.csv','.png'): prefix=prefix.with_suffix('')
            paths=[Path(str(prefix)+ext) for ext in ('.json','.csv','.png')]
            if self.path and any(p.resolve()==self.path.resolve() for p in paths):
                raise ValueError('输出名称会覆盖输入文件，请换一个名称。')
            if not overwrite and any(p.exists() for p in paths):
                if not messagebox.askyesno('覆盖结果','该名称的部分结果文件已存在，是否覆盖？'): return
            prefix.parent.mkdir(parents=True,exist_ok=True)
            r=self.result
            payload=dict(frequency_unit='THz',power_conversion='amplitude**2, power unchanged, dB:10**(dB/10)',
                         input=self.snapshot,result=r.to_dict())
            paths[0].write_text(json.dumps(clean_json(payload),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
            with paths[1].open('w',encoding='utf-8-sig',newline='') as fh:
                writer=csv.writer(fh); writer.writerow(['frequency_THz','power','fitted_power','baseline','residual'])
                for i in range(len(r.f)):
                    fit=float(r.fitted[i]) if len(r.fitted) else ''
                    base=float(r.baseline[i]) if len(r.baseline) else ''
                    writer.writerow([r.f[i],r.power[i],fit,base,float(r.power[i]-fit) if fit!='' else ''])
            self.fig.savefig(paths[2],dpi=180)
            self.status.set(f'已保存：{prefix}.json / .csv / .png')
        except Exception as exc: self.error(exc)

    def error(self,exc):
        self.status.set('操作未完成：'+str(exc)); messagebox.showerror('请检查设置',str(exc))


def main():
    parser=argparse.ArgumentParser(description='图形界面谐振 Q 计算器')
    parser.add_argument('file',nargs='?'); args=parser.parse_args()
    root=tk.Tk(); app=QCalculator(root)
    if args.file: root.after(100,lambda:app.load_file(args.file))
    root.mainloop()


if __name__=='__main__': main()

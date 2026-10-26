"""Generate synthetic examples (not measured spectra) and fit all five models."""
from pathlib import Path
import numpy as np
from scipy.special import voigt_profile
from qfit import fit_resonance, save_result

out = Path('results/qfit_demo')
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(2026)
f = np.linspace(.27, .39, 1201)
f0, width = .33, .006
eps = 2 * (f-f0) / width
lorentz = 1 / (1+eps**2)
gaussian = np.exp(-4*np.log(2)*((f-f0)/width)**2)
sigma, gamma = .004/np.sqrt(8*np.log(2)), .003/2
v = voigt_profile(f-f0, sigma, gamma) / voigt_profile(0, sigma, gamma)
curves = {
    'lorentzian': .12 + .75*lorentz,
    'fano': .4 + .25*(((1.5+eps)**2/(1+eps**2))-1)/3.25,
    'gaussian': .12 + .75*gaussian,
    'voigt': .12 + .75*v,
    'pseudo_voigt': .12 + .75*(.35*lorentz+.65*gaussian),
}
for model, power in curves.items():
    power = power + rng.normal(0, .0005, len(f))
    source = out / (model + '_synthetic.txt')
    np.savetxt(source, np.column_stack([f, np.sqrt(power)]), header='frequency_THz amplitude_SYNTHETIC')
    fit = fit_resonance(f, np.sqrt(power), model=model, feature='peak')
    save_result(fit, out / model)
    print(model, 'pole_Q=', fit.q_resonance, 'spectral_Q=', fit.q_spectral)

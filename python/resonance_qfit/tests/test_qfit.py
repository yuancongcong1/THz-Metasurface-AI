import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.special import voigt_profile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qfit import fit_resonance, load_spectrum, save_result, to_power
from batch_qfit import main as batch_main


class ResonanceFits(unittest.TestCase):
    def setUp(self):
        self.f = np.linspace(.27, .39, 601)
        self.f0, self.gamma = .33, .006
        self.eps = 2 * (self.f - self.f0) / self.gamma

    def test_power_conversion(self):
        np.testing.assert_allclose(to_power([.1, .5, 1], 'amplitude'), [.01, .25, 1])
        np.testing.assert_allclose(to_power([-20, -3.01029995664, 0], 'db'), [.01, .5, 1], rtol=1e-10)
        np.testing.assert_allclose(to_power([.3 + .4j], 'complex'), [.25])
        with self.assertRaises(ValueError): to_power([-.1], 'amplitude')

    def test_lorentzian_and_units(self):
        p = .15 + .8 / (1 + self.eps ** 2) + .1 * (self.f - .33)
        results = []
        for factor, unit in [(1, 'THz'), (1e3, 'GHz'), (1e12, 'Hz')]:
            r = fit_resonance(self.f * factor, np.sqrt(p), model='lorentzian', frequency_unit=unit, n_starts=6)
            results.append(r)
            self.assertAlmostEqual(r.q_resonance, 55, places=3)
            self.assertAlmostEqual(r.q_spectral, 55, places=2)
            self.assertLess(r.diagnostics['rmse_power'], 1e-7)
        self.assertAlmostEqual(results[0].q_resonance, results[2].q_resonance, places=4)

    def test_fano_pole_not_feature_q(self):
        q = 1.5
        p = .4 + .25 * (((q + self.eps) ** 2 / (1 + self.eps ** 2)) - 1) / (1 + q*q)
        r = fit_resonance(self.f, p, kind='power', model='fano', feature='peak', n_starts=8)
        self.assertAlmostEqual(r.f0, .33, places=6)
        self.assertAlmostEqual(r.q_resonance, 55, places=3)
        self.assertAlmostEqual(r.parameters['fano_q_canonical_positive_scale'], q, places=3)
        self.assertAlmostEqual(r.feature_frequency, .332, places=5)
        self.assertGreater(abs(r.q_spectral - r.q_resonance), 5)

    def test_negative_fano_q_and_dip(self):
        q = -.8
        p = .5 + .2 * (((q + self.eps)**2 / (1 + self.eps**2)) - 1) / (1 + q*q)
        r = fit_resonance(self.f, p, kind='power', model='fano', feature='dip', n_starts=8)
        self.assertAlmostEqual(r.q_resonance, 55, places=3)
        self.assertAlmostEqual(r.parameters['fano_q_canonical_positive_scale'], q, places=3)
        self.assertIsNone(r.q_3db)
        self.assertGreater(r.q_spectral, 0)

    def test_lorentzian_dip(self):
        p = .9 - .7 / (1 + self.eps ** 2)
        r = fit_resonance(self.f, p, kind='power', model='lorentzian', feature='dip', n_starts=6)
        self.assertAlmostEqual(r.q_resonance, 55, places=3)
        self.assertAlmostEqual(r.q_spectral, 55, places=2)

    def test_absolute_half_power_is_not_half_amplitude(self):
        power = .8 / (1+self.eps**2)
        r = fit_resonance(self.f, np.sqrt(power), model='lorentzian', n_starts=6)
        self.assertAlmostEqual(r.q_3db, 55, places=2)

    def test_gaussian_q_is_spectral_only(self):
        p = .1 + .8 * np.exp(-4 * np.log(2) * ((self.f - .33)/.006)**2)
        r = fit_resonance(self.f, p, kind='power', model='gaussian', n_starts=6)
        self.assertIsNone(r.q_resonance)
        self.assertAlmostEqual(r.q_spectral, 55, places=2)

    def test_pseudo_voigt(self):
        g = np.exp(-4 * np.log(2) * ((self.f - .33)/.006)**2)
        p = .1 + .8 * (.35 / (1 + self.eps**2) + .65*g)
        r = fit_resonance(self.f, p, kind='power', model='pseudo_voigt', n_starts=8)
        self.assertAlmostEqual(r.q_spectral, 55, places=2)
        self.assertAlmostEqual(r.parameters['lorentz_fraction'], .35, places=3)

    def test_voigt(self):
        sigma, gamma = .004 / np.sqrt(8 * np.log(2)), .003 / 2
        p = .12 + .8 * voigt_profile(self.f-.33, sigma, gamma) / voigt_profile(0, sigma, gamma)
        r = fit_resonance(self.f, p, kind='power', model='voigt', n_starts=8)
        self.assertAlmostEqual(r.parameters['gaussian_fwhm'], .004, places=6)
        self.assertAlmostEqual(r.parameters['lorentzian_fwhm'], .003, places=6)
        self.assertIsNone(r.q_resonance)
        self.assertAlmostEqual(r.q_spectral, r.f0 / r.linewidth, places=3)

    def test_noisy_auto_selection(self):
        rng = np.random.default_rng(21)
        p = .15 + .7 / (1+self.eps**2) + rng.normal(0, .001, len(self.f))
        r = fit_resonance(self.f, p, kind='power', model='auto', n_starts=6)
        self.assertAlmostEqual(r.f0, .33, places=4)
        self.assertLess(abs(r.q_spectral/55 - 1), .02)
        self.assertEqual(len(r.candidates), 5)

    def test_undersampling_and_truncation_flags(self):
        f = np.linspace(.27, .39, 25)
        p = .1 + .8 / (1 + (2*(f-.33)/.006)**2)
        r = fit_resonance(f, p, kind='power', model='lorentzian', n_starts=6)
        self.assertTrue(any('Undersampled' in w for w in r.warnings))
        f = np.linspace(.33, .37, 101)
        p = .1 + .8 / (1 + (2*(f-.33)/.006)**2)
        r = fit_resonance(f, p, kind='power', model='lorentzian', n_starts=6)
        self.assertIsNone(r.feature_fwhm)

    def test_input_rejection_and_sorting(self):
        p = .1 + .8 / (1 + self.eps**2)
        for badf, badp in [(self.f[:5], p[:5]), (self.f, np.ones(len(self.f))),
                            (self.f, .2+.1*self.f),
                            (np.r_[self.f[:-1], self.f[-2]], p), (self.f, p * np.nan)]:
            with self.assertRaises(ValueError): fit_resonance(badf, badp, kind='power')
        r = fit_resonance(self.f[::-1], p[::-1], kind='power', model='lorentzian', n_starts=4)
        self.assertAlmostEqual(r.q_resonance, 55, places=3)

    def test_cst_import_and_exports(self):
        p = .1 + .8 / (1 + self.eps**2)
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'spectrum.txt'
            with path.open('w') as handle:
                handle.write('# CST example: Frequency / THz; magnitude\n----------\n')
                np.savetxt(handle, np.column_stack([self.f, np.sqrt(p)]))
            f, a = load_spectrum(path)
            r = fit_resonance(f, a, model='lorentzian', n_starts=4)
            prefix = Path(td) / 'fit'
            save_result(r, prefix)
            loaded = json.loads(prefix.with_suffix('.json').read_text(encoding='utf-8'))
            self.assertAlmostEqual(loaded['q_resonance'], 55, places=3)
            self.assertTrue(prefix.with_suffix('.png').is_file())
            self.assertTrue(prefix.with_suffix('.csv').is_file())

    def test_batch_records_failures(self):
        with tempfile.TemporaryDirectory() as td:
            folder = Path(td) / 'input'; folder.mkdir()
            out = Path(td) / 'fits'
            p = .1+.8/(1+self.eps**2)
            np.savetxt(folder / 'good.txt', np.column_stack([self.f, np.sqrt(p)]))
            np.savetxt(folder / 'bad.txt', [[.33,.4],[.33,.5]])
            code = batch_main([str(folder), '--window', '.27', '.39', '--model',
                               'lorentzian', '--output-dir', str(out), '--no-plots'])
            self.assertEqual(code, 1)
            self.assertTrue((out / 'good.txt_fit.json').is_file())
            self.assertIn('Duplicate frequencies', (out / 'summary.csv').read_text(encoding='utf-8-sig'))


if __name__ == '__main__':
    unittest.main()

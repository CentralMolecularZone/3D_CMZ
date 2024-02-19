from __future__ import division
import numpy as np
import pyspeckit
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt

# Initial guesses for pyspeckit Gaussian fitting
amp_guess = {1: '0.12', 2: '1', 3: '0.05', 4: '0.35', 5: '0.28', 6: '1', 7: '1',
             8: '1.3', 9: '0.9', 10: '1.25', 11: '1.2', 12: '1.2', 13: '1.9',
             14: '2.4', 15: '0.25', 16: '0.42', 17: '1.3', 18: '1.2', 19: '2.4',
             20: '0.9', 21: '0.25', 22: '0.25'}

v_cen_guess = {1: '19', 2: '39', 3: '16', 4: '-56', 5: '-25', 6: '58', 7: '84',
               8: '15', 9: '54', 10: '83', 11: '50', 12: '48', 13: '50',
               14: '62', 15: '35', 16: '52', 17: '49', 18: '29', 19: '53',
               20: '20', 21: '39', 22: '-2'}

fwhm_guess = {1: '22', 2: '34', 3: '29', 4: '22', 5: '24', 6: '26', 7: '20',
              8: '27', 9: '22', 10: '20', 11: '22', 12: '25', 13: '15', 14: '29',
              15: '32', 16: '23', 17: '23', 18: '25', 19: '28', 20: '25',
              21: '10', 22: '30'}

amp_guess2 = {1: '0.1', 6: '1', 7: '0.3', 15: '0.25', 21: '0.25'}

v_cen_guess2 = {1: '-20', 6: '28', 7: '74', 15: '35', 21: '8'}

fwhm_guess2 = {1: '30', 6: '35', 7: '38', 15: '32', 21: '30'}

# Loop through all leaves to fit HNCO spectra
for i in range(1,23):
    if i in amp_guess2:
        guess = [amp_guess[i], v_cen_guess[i], fwhm_guess[i]] + [amp_guess2[i],
                v_cen_guess2[i], fwhm_guess2[i]]
        sp = pyspeckit.Spectrum('./../Leaf_cubes_MALT90/meanspec/'+str(i)+'_HNCO_cube_meanspec.fits')
        sp.plotter(linestyle='--')
        sp.specfit(fittype='gaussian', Interactive=False, color='red',
                   guesses=guess, annotate=False)
        amp  = np.around(sp.specfit.parinfo['AMPLITUDE0'].value, decimals=2)
        vel  = np.around(sp.specfit.parinfo['SHIFT0'].value, decimals=1)
        sig  = np.around(sp.specfit.parinfo['WIDTH0'].value, decimals=1)
        fwhm = np.around(sp.specfit.parinfo['WIDTH0'].value * np.sqrt(8*np.log(2)),
                         decimals=1)
        amp2  = np.around(sp.specfit.parinfo['AMPLITUDE1'].value, decimals=2)
        vel2  = np.around(sp.specfit.parinfo['SHIFT1'].value, decimals=1)
        sig2  = np.around(sp.specfit.parinfo['WIDTH1'].value, decimals=1)
        fwhm2 = np.around(sp.specfit.parinfo['WIDTH1'].value * np.sqrt(8*np.log(2)),
                          decimals=1)
        plt.annotate("Peak   = [" + str(amp) +", "+ str(amp2) + "] K",
                     xy=(0.05, 0.95), xycoords='axes fraction')
        plt.annotate("$v_{\mathrm{cen}}$     = [" + str(vel) +", "+ str(vel2) + "] km/s",
                     xy=(0.05, 0.9), xycoords='axes fraction')
        plt.annotate("FWHM = [" + str(fwhm) +", "+ str(fwhm2) + "] km/s",
                     xy=(0.05, 0.85), xycoords='axes fraction')
        plt.annotate("$\sigma$         = [" + str(sig) +", "+ str(sig2) + "] km/s",
                     xy=(0.05, 0.8), xycoords='axes fraction')
        plt.title("Fitted averaged HNCO spectrum for structure "+str(i))
        plt.tight_layout()
        sp.plotter.savefig('./../Figs/HNCO_fits/'+str(i)+'_HNCO_fit.pdf')
        plt.close()
    else:
        guess = [amp_guess[i], v_cen_guess[i], fwhm_guess[i]]
        sp = pyspeckit.Spectrum('./../Leaf_cubes_MALT90/meanspec/'+str(i)+'_HNCO_cube_meanspec.fits')
        sp.plotter(linestyle='--')
        sp.specfit(fittype='gaussian', Interactive=False, color='red',
                   guesses=guess, annotate=False)
        amp  = np.around(sp.specfit.parinfo['AMPLITUDE0'].value, decimals=2)
        vel  = np.around(sp.specfit.parinfo['SHIFT0'].value, decimals=1)
        sig  = np.around(sp.specfit.parinfo['WIDTH0'].value, decimals=1)
        fwhm = np.around(sp.specfit.parinfo['WIDTH0'].value * np.sqrt(8*np.log(2)), decimals=1)
        plt.annotate("Peak   = " + str(amp) + " K", xy=(0.05, 0.95), xycoords='axes fraction')
        plt.annotate("$v_{\mathrm{cen}}$     = " + str(vel) + " km/s", xy=(0.05, 0.9), xycoords='axes fraction')
        plt.annotate("FWHM = " + str(fwhm) + " km/s", xy=(0.05, 0.85), xycoords='axes fraction')
        plt.annotate("$\sigma$         = " + str(sig) + " km/s", xy=(0.05, 0.8), xycoords='axes fraction')
        plt.title("Fitted averaged HNCO spectrum for structure "+str(i))
        plt.tight_layout()
        sp.plotter.savefig('./../Figs/HNCO_fits/'+str(i)+'_HNCO_fit.pdf')
        plt.close()

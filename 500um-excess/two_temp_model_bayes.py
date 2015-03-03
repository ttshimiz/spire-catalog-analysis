"""
Module defining functions to fit a two temperature
modified blackbody (greybody) model.
The fitting routine uses a Monte Carlo Markov Chain to sample the posterior
distribution. 2.5, 16, 50, 84, and 97.5 %tiles are determined for each
model parameter.
Currently only N_cold, T_cold, Beta_cold, and N_warm are free parameters with
T_warm and Beta_warm fixed at 60 and 2.0 respectively.
The model calculates simulated photometry using the transmission filters for
the Herschel and WISE instruments as well as the given redshift of the source.

Priors
------
N_cold: Uniform between -inf and inf
T_cold: Uniform between 0 and 50 K
Beta_cold: Uniform between 0 and 5
N_warm: Uniform between -inf and inf

The likelihood function is a gaussian.
"""

import numpy as np
import emcee
import astropy.constants as const
import astropy.units as u
from astropy.units.quantity import Quantity

# GLOBAL CONSTANTS
planck_h = const.h.cgs.value
k_b = const.k_B.cgs.value
c = const.c.cgs.value
pc2cm = const.pc.cgs.value
Msun2gram = const.M_sun.cgs.value
c_micron = const.c.to('micron / s').value

# FILTERS
w4_filter = np.loadtxt('Filters/WISE/RSR-W4.EE.txt')
h70_filter = np.loadtxt('Filters/Herschel/PACS70.txt')
h70_ind = (h70_filter[:, 0] > 50) & (h70_filter[:, 0] < 105)
h70_filter = h70_filter[h70_ind]
h160_filter = np.loadtxt('Filters/Herschel/PACS160.txt')
h250_filter = np.loadtxt('Filters/Herschel/PSW.txt')
sort250 = np.argsort(h250_filter[:, 0])
h250_filter = h250_filter[sort250]
h250_ind = (h250_filter[:, 0] > 180) & (h250_filter[:, 0] < 310)
h250_filter = h250_filter[h250_ind]
h350_filter = np.loadtxt('Filters/Herschel/PMW.txt')
sort350 = np.argsort(h350_filter[:, 0])
h350_filter = h350_filter[sort250]
h350_ind = (h350_filter[:, 0] > 250) & (h350_filter[:, 0] < 450)
h350_filter = h350_filter[h350_ind]
h500_filter = np.loadtxt('Filters/Herschel/PLW.txt')
sort500 = np.argsort(h500_filter[:, 0])
h500_filter = h500_filter[sort250]
h500_ind = (h500_filter[:, 0] > 360) & (h500_filter[:, 0] < 620)
h500_filter = h500_filter[h500_ind]


# Class to contain SED information
class SED(object):

    def __init__(self, waves, flux, error='None', z=0.0):
        self.set_wavelength(waves)
        self.set_flux(flux)
        if (error == 'None'):
            self._error = np.ones(len(self.flux)) * self.flux.unit
        else:
            self.set_error(error)
        self.set_redshift(z)

    def set_wavelength(self, waves):
        if isinstance(waves, Quantity):
            self.waves = waves
        else:
            raise TypeError('Wavelengths must have units attached to it!')

    def set_flux(self, flux):
        if isinstance(flux, Quantity):
            self.flux = flux
        else:
            raise TypeError('Fluxes must have units attached to it!')

    def set_error(self, error):
        if isinstance(error, Quantity):
            if (error.unit == self.flux.unit):
                self.error = error
            else:
                raise ValueError('Units on error must be same as'
                                 'units on flux!')
        else:
            raise TypeError('Errors must have units attached to it!')

    def set_redshift(self, z):
        self.redshift = z

    def waves_to(self, unit):
        if isinstance(unit, u.core.Unit):
            self.waves = self.waves.to(unit)
        else:
            raise TypeError('"unit" must be an Astropy Unit instance')

    def flux_to(self, unit):
        if isinstance(unit, u.core.Unit):
            self.flux = self.flux.to(unit)
            self.error = self.error.to(unit)
        else:
            raise TypeError('"unit" must be an Astropy Unit instance')


# Routine to convolve a transmission curve with an SED to get a
# monochromatic flux density
def convolve_with_trans(sed_waves, sed_flux, filt_waves, filt_trans):
    """
    Function to integrate the SED over a bandpass with a defined transmission
    curve

    Input Parameters
    ----------------
        sed_waves   = wavelengths at which the SED is defined
        sed_flux    = flux densities of the SED
        filt_waves  = wavelengths at which the filter is defined
        filt_trans  = transmission of the filter

    Output
    ------
        fnu = integrated flux density

    """

    # Interpolate the transmission to the same wavelengths as the SED
    # Use a value of 0 for where the filter isn't defined
    if ~all(sed_waves == filt_waves):
        interp_trans = np.interp(sed_waves, filt_waves, filt_trans,
                                 left=0, right=0)
    else:
        interp_trans = filt_trans

    # Integrate over all wavelengths (but use the frequencies
    # since the SED is in Jy)
    integ_num = np.trapz(c_micron/sed_waves, interp_trans*sed_flux)
    integ_denom = np.trapz(c_micron/sed_waves, interp_trans)

    return integ_num/integ_denom


# Routine for generating the model flux densities for a modified blackbody
def greybody(nu, norm, tdust, beta):
    """
    Single temperature greybody component of the model to
    describe the cold dust
    Equation: S_grey = mdust * kappa0 * (nu / (c/lamdda_norm))^beta *
                       B(nu, tdust) / lumD^2

    Input Parameters
    ----------------
       nu           = frequencies at which to evaluate the greybody in Hz
       norm         = Normalization of the greybody
       tdust        = dust temperature in Kelvin
       beta         = emissivity describing the dependence of the absorption
                      coefficient on frequency

    Output
    ------
       snu = array of flux densities in Jy with same shape as waves

    """
    aa = nu**beta
    b_nu_A = (2 * planck_h * nu**3 / c**2)
    b_nu_B = 1 / (np.exp(planck_h * nu / k_b / tdust) - 1)
    b_nu = b_nu_A * b_nu_B
    snu = 10**norm * aa * np.pi * b_nu

    # Convert to Jy
    snu = 10**(np.log10(snu) + 23)

    return snu


# Combine the two components together
def twotemp_model(waves, norm_cold, tdust_cold, beta_cold,
                  norm_warm, tdust_warm=60., beta_warm=2.):
    mbb_cold = greybody(c_micron/waves, norm_cold, tdust_cold, beta_cold)
    mbb_warm = greybody(c_micron/waves, norm_warm, tdust_warm, beta_warm)

    return mbb_cold+mbb_warm


# Calculate monochromatic flux densities of the model using the filter
# transmission curves and redshift of the source
def calc_model(waves, norm_cold, tdust_cold, beta_cold,
               norm_warm, tdust_warm, beta_warm, redshift):
    """
    For each observed wavelength we need to determine which filter to use.
    22. micron will be assumed to be WISE filters.
    70 - 500 micron will be assumed to be Herschel.
    Then the wavelengths in the filter need to be redshifted to get the
    rest frame wavelengths.
    The model will be calculated at these rest frame wavelengths,
    then redshifted back to observed frame
    where it will be convolved with the transmission curve.
    """

    # Array to store the calculated model fluxes
    model_fluxes = np.zeros(np.shape(waves))

    for i, w in enumerate(waves):
        if w == 22.:
            filter = w4_filter
        elif w == 70.:
            filter = h70_filter
        elif w == 160.:
            filter = h160_filter
        elif w == 250.:
            filter = h250_filter
        elif w == 350.:
            filter = h350_filter
        elif w == 500.:
            filter = h500_filter

        filter_waves = filter[:, 0]
        filter_trans = filter[:, 1]

        # Transform to rest frame
        rest_waves = filter_waves / (1 + redshift)

        # Calculate the model at the rest frame wavelengths
        rest_total_model = twotemp_model(rest_waves, norm_cold, tdust_cold,
                                         beta_cold, norm_warm, tdust_warm,
                                         beta_warm)

        # Convolve the SED with the transmission curve
        mono_flux = convolve_with_trans(filter_waves, rest_total_model,
                                        filter_waves, filter_trans)

        model_fluxes[i] = mono_flux

    return model_fluxes


# Log-likelihood
def log_like(theta, x, y, sigma, fix_beta_cold, tdust_warm, beta_warm, zz):
    if ~fix_beta_cold:
        norm_cold, tdust_cold, beta_cold, norm_warm = theta
    else:
        norm_cold, tdust_cold, norm_warm = theta
        beta_cold = fix_beta_cold
    model = calc_model(x, norm_cold, tdust_cold, beta_cold,
                       norm_warm, tdust_warm, beta_warm, zz)

    return -0.5*(np.sum((y - model)**2 / sigma**2 +
                 np.log(2 * np.pi * sigma**2)))


# Log prior distribution, uniform priors with limits for beta_cold, tdust_cold,
# norm_cold, and norm_warm
def log_prior(theta, fix_beta_cold):

    if ~fix_beta_cold:
        norm_cold, tdust_cold, beta_cold, norm_warm = theta
        if 0 < tdust_cold < 50 and 0 < beta_cold < 5:
            return 0
    else:
        norm_cold, tdust_cold, norm_warm = theta
        if 0 < tdust_cold < 50:
            return 0
    return -np.inf


# Log posterior = log-likelihood + log-prior
def log_post(theta, x, y, yerr, fix_beta_cold, tdust_warm, beta_warm, zz):
    lp = log_prior(theta, fix_beta_cold)
    llike = log_like(theta, x, y, yerr,
                     fix_beta_cold, tdust_warm, beta_warm, zz)
    if not np.isfinite(lp) or not np.isfinite(llike):
        return -np.inf
    return lp + llike


# Main fitting routine
def fit_two_temp_bayes(sed, fix_beta_cold=False, tdust_warm=60,
                       beta_warm=2, nwalkers=50, nburn=200, nsteps=1000):

    # Convert fluxes to Jy
    sed.flux_to(u.Jy)

    # Convert wavelengths to micron
    sed.waves_to(u.micron)

    # x values are wavelengths
    # y values are fluxes
    # y errors are photometric errors
    x = sed.waves.value
    y = sed.flux.value
    yerr = sed.error.value
    zz = sed.redshift

    # Initial guesses for parameters
    tdust_warm_guess = tdust_warm
    beta_warm_guess = beta_warm
    norm_warm_guess = np.log10(y[0]/greybody(c_micron/22., 0.,
                               tdust_warm_guess, beta_warm_guess))
    tdust_cold_guess = 25.
    if ~fix_beta_cold:
        beta_cold_guess = 2.0
        ndims = 4
    else:
        beta_cold_guess = fix_beta_cold
        ndims = 3
    norm_cold_guess = np.log10(y[2]/greybody(c_micron/160., 0.,
                               tdust_cold_guess, beta_cold_guess))

    # Need guesses for each walker
    # Randomly perturb from initial guesses
    if ~fix_beta_cold:
        guess = np.array([[norm_cold_guess, tdust_cold_guess, beta_cold_guess,
                           norm_warm_guess] + 1e-4*np.random.randn(ndims)
                          for k in range(nwalkers)])
    else:
        guess = np.array([[norm_cold_guess, tdust_cold_guess, norm_warm_guess]
                          + 1e-4*np.random.randn(ndims)
                          for k in range(nwalkers)])

    # Setup MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndims, log_post,
                                    args=(x, y, yerr, fix_beta_cold,
                                          tdust_warm, beta_warm, zz),
                                    threads=4)
    # Run the MCMC
    sampler.run_mcmc(guess, nsteps)

    # Extract the chains and flatten into one long chain
    samples_all = sampler.chain[:, :, :].reshape(-1, ndims)
    samples_noburn = sampler.chain[:, nburn:, :].reshape(-1, ndims)

    # Close the Pool
    sampler.pool.close()

    # Calculate the 2.5, 16, 50, 84, and 97.5 percentiles
    norm_cold = np.percentile(samples_noburn[:, 0], [2.5, 16, 50, 84, 97.5])
    tdust_cold = np.percentile(samples_noburn[:, 1], [2.5, 16, 50, 84, 97.5])
    beta_cold = np.percentile(samples_noburn[:, 2], [2.5, 16, 50, 84, 97.5])
    norm_warm = np.percentile(samples_noburn[:, 3], [2.5, 16, 50, 84, 97.5])

    # Return the percentiles and the samples
    result = {'samples_all': samples_all,
              'samples_noburn': samples_noburn,
              'tdust_warm': tdust_warm_guess,
              'beta_warm': beta_warm_guess,
              'norm_warm': norm_warm,
              'norm_cold': norm_cold,
              'tdust_cold': tdust_cold,
              'beta_cold': beta_cold,
              'nwalkers': nwalkers,
              'nburn': nburn,
              'nsteps': nsteps}

    return result

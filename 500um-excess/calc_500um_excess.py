# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 17:22:48 2015

Script to run the two temperature fitting on the BAT AGN that have detections
at all the wavelengths from 22 - 500 microns. However only the 22-350 micron
points will be fit.
After fitting, the 500 micron excess will be measured by comparing the best fit
model 500 micron point with the observed 500 micron point.
The excess will be measured as a percentage in the following way:
    excess500 = [F_obs(500) - F_model(500)] / F_model(500)

The excess values along with fit results will be saved in a text file.

@author: ttshimiz
"""

import numpy as np
import pandas as pd
import two_temp_model_bayes as ttmb
import astropy.units as u
import triangle
import matplotlib.pyplot as plt
import seaborn
import pickle

seaborn.set()

# Directory where all of the data is stored
data_dir = '/Users/ttshimiz/Research/Thesis/bat-data/'


# Functions for analysis
def excess(model, observed):
    return (observed - model) / model


def calc_excess(trace, obs500, zz=0, fix_beta_cold=False, fix_tdust_warm=60.,
                beta_warm=2.0):
    """
    Function that calculates the 500 micron excess for every modeled SED.

    INPUTS
    ------
    trace = Array of parameter values from the MCMC fitting routine that
            defines the posterior distribution
    obs500 = The observed 500 micron flux in Jy
    zz = redshift of the source (Default to 0)
    fix_beta_cold = False or a scalar value. If False then the emissivity for
                    the cold component was allowed to vary. Scalar value gives
                    the emissivity if it was fixed. (Default False)
    tdust_warm = The fixed temperature of the warm component. (Default 60 K)
    beta_warm = The fixed emissivity of the warm component. (Default 2)

    OUTPUTS
    -------
    excess500 = Dictionary containing the 500 micron excess for each model as
                well as the calculated mean, 16th, 50th, and 84th percentile.
    """

    # Calculate the modeled 500 micron flux
    norm_cold_model = trace[:, 0]
    tdust_cold_model = trace[:, 1]
    if (fix_beta_cold is None):
        beta_cold_model = trace[:, 2]
        norm_warm_model = trace[:, 3]
        tdust_warm_model = np.ones(len(norm_cold_model))*fix_tdust_warm
    elif (fix_tdust_warm is None):
        beta_cold_model = np.ones(len(norm_cold_model))*fix_beta_cold
        norm_warm_model = trace[:, 2]
        tdust_warm_model = trace[:, 3]
    else:
    	beta_cold_model = np.ones(len(norm_cold_model))*fix_beta_cold
        tdust_warm_model = np.ones(len(norm_cold_model))*fix_tdust_warm
        norm_warm_model = trace[:, 2]

    model500 = np.zeros(len(norm_cold_model))
    for i in range(len(norm_cold_model)):
        model500[i] = ttmb.calc_model([500.], norm_cold_model[i],
                                      tdust_cold_model[i],
                                      beta_cold_model[i], norm_warm_model[i],
                                      tdust_warm_model[i], beta_warm, zz)

    e500 = excess(model500, obs500)
    p2_5, p16, p50, p84, p97_5 = np.percentile(e500, q=[2.5, 16, 50, 84, 97.5])
    excess500 = {'all': e500,
                 'mean': np.mean(e500),
                 'p2_5': p2_5,
                 'p16': p16,
                 'p50': p50,
                 'p84': p84,
                 'p97_5': p97_5}

    return excess500


def run_all_sources(fix_beta_cold=None, fix_tdust_warm=60., beta_warm=2.0,
                    nwalkers=50, nburn=200, nsteps=1000,
                    results_dir='results/betaFree/',
                    fsuffix='_betafree_twotemp_500excess'):

    # Upload the Herschel and WISE data
    herschel_data = pd.read_csv(data_dir+'bat_herschel.csv', index_col=0)
    wise_data = pd.read_csv(data_dir+'bat_wise.csv', index_col=0)

    # Upload the general information for the BAT sources
    gen_data = pd.read_csv(data_dir+'bat_info.csv', index_col=0)

    # We will use the 22, 70, 160, 250, 350, and 500 micron data
    # We also need the redshift as well
    fluxes = pd.DataFrame({'W4': wise_data.W4,
                           'H70': herschel_data.H_70,
                           'H160': herschel_data.H_160,
                           'H250': herschel_data.H_250,
                           'H350': herschel_data.H_350,
                           'H500': herschel_data.H_500})

    flux_err = pd.DataFrame({'W4_err': wise_data.W4_err,
                             'H70_err': herschel_data.H_70_err,
                             'H160_err': herschel_data.H_160_err,
                             'H250_err': herschel_data.H_250_err,
                             'H350_err': herschel_data.H_350_err,
                             'H500_err': herschel_data.H_500_err})

    redshifts = gen_data.Redshift

    # Find all of the sources that are detected at every wavelength
    detect_all = (fluxes > 0).all(axis=1)

    # Iterate over all of the sources
    waves = np.array([22., 70., 160., 250., 350.])*u.micron

    for idx in fluxes[detect_all].index[131:]:
        print 'Running analysis on ...' + idx
        flux_src = fluxes.loc[idx][['W4', 'H70', 'H160', 'H250', 'H350']]
        flux_err_src = flux_err.loc[idx][['W4_err', 'H70_err', 'H160_err',
                                          'H250_err', 'H350_err']]
        zz_src = redshifts.loc[idx]

        sed = ttmb.SED(waves,
                       np.array(flux_src)*u.Jy,
                       np.array(flux_err_src)*u.Jy,
                       z=zz_src)
        print '\tFitting model...'
        fit_results = ttmb.fit_two_temp_bayes(sed, fix_beta_cold=fix_beta_cold,
                                              fix_tdust_warm=fix_tdust_warm,
                                              beta_warm=beta_warm,
                                              nwalkers=nwalkers,
                                              nburn=nburn,
                                              nsteps=nsteps)

        # Calculate the 500 micron excess
        flux500 = fluxes.loc[idx]['H500']
        print '\tCalculating 500 micron excess...'
        excess500 = calc_excess(fit_results['samples_noburn'], flux500, zz_src,
                                fix_beta_cold=fix_beta_cold,
                                fix_tdust_warm=fix_tdust_warm)
        print '\tPlotting and saving results...'

        if (fix_beta_cold is None):
            labels = [r'$N_{\mathrm{cold}}$',
                      r'$T_{\mathrm{cold}}$',
                      r'$\beta_{\mathrm{cold}}$',
                      r'$N_{\mathrm{warm}}$']
        elif (fix_tdust_warm is None):
            labels = [r'$N_{\mathrm{cold}}$',
                      r'$T_{\mathrm{cold}}$',
                      r'$N_{\mathrm{warm}}$',
                      r'$T_{\mathrm{warm}}$']
        else:
            labels = [r'$N_{\mathrm{cold}}$',
                      r'$T_{\mathrm{cold}}$',
                      r'$N_{\mathrm{warm}}$']

        # Plot a corner triangle figure to look at relationships
        # between parameters
        fig_tri = triangle.corner(fit_results['samples_noburn'],
                                  labels=labels,
                                  quantiles=[0.16, 0.5, 0.84],
                                  verbose=False)

        # Save the triangle plot
        fig_tri.savefig(results_dir + 'triangle_plots/' + idx +
                        fsuffix + '_triangle.png',
                        bbox_inches='tight')
        plt.close(fig_tri)

        # Plot the data along with the median SED. Shade in the 95% confidence
        # interval
        # Wavelengths to calculate the model SEDs at
        model_wave = np.arange(1, 1000, 1)

        # Randomly choose 1000 modeled parameters
        p_rand = np.random.randint(low=0, high=((nsteps-nburn)*nwalkers),
                                   size=1000)
        model_ncold = fit_results['samples_noburn'][p_rand, 0][:, None]
        model_tcold = fit_results['samples_noburn'][p_rand, 1][:, None]

        if (fix_beta_cold is None):
            model_betacold = fit_results['samples_noburn'][p_rand, 2][:, None]
            model_nwarm = fit_results['samples_noburn'][p_rand, 3][:, None]
            model_twarm = (np.ones(len(model_ncold))*fix_tdust_warm)[:, None]

        elif (fix_tdust_warm is None):
            model_betacold = (np.ones(len(model_ncold))*fix_beta_cold)[:, None]
            model_nwarm = fit_results['samples_noburn'][p_rand, 2][:, None]
            model_twarm = fit_results['samples_noburn'][p_rand, 3][:, None]

        else:
            model_betacold = (np.ones(len(model_ncold))*fix_beta_cold)[:, None]
            model_twarm = (np.ones(len(model_ncold))*fix_tdust_warm)[:, None]
            model_nwarm = fit_results['samples_noburn'][p_rand, 2][:, None]

        model_betawarm = (np.ones(len(model_ncold))*beta_warm)[:, None]
        model_fluxes = ttmb.twotemp_model(model_wave/(1+zz_src),
                                          model_ncold,
                                          model_tcold,
                                          model_betacold,
                                          model_nwarm,
                                          model_twarm,
                                          model_betawarm)

        # 2.5%, 50%, and 97.5% percentile to form median and
        # 95% spread in model SEDs
        model_2_5, model_50, model_97_5 = np.percentile(model_fluxes,
                                                        [2.5, 50., 97.5],
                                                        axis=0)

        # Plotting SEDs
        cp = seaborn.color_palette()
        fig_sed = plt.figure()
        ax = fig_sed.add_subplot(1, 1, 1)
        ax.loglog(model_wave, model_50, color=cp[0], lw=2)
        ax.fill_between(model_wave, model_2_5, model_97_5, color=cp[0],
                        alpha=0.3)
        ax.errorbar(waves.value, flux_src, yerr=flux_err_src, fmt='o',
                    color='k', ms=8, ls='None')
        ax.errorbar(500., flux500, yerr=flux_err.loc[idx]['H500_err'],
                    fmt='o', color=cp[2], ms=8, ls='None')
        ax.set_ylim(np.min([flux_src.min(), flux500])*10**(-0.5),
                    np.max(model_97_5)*10**(0.5))
        ax.set_xlim(10, 1000)
        ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        ax.set_ylabel(r'$F_{\nu}$ [Jy]')
        ax.text(0.05, 0.95, idx, ha='left', va='center',
                transform=ax.transAxes)

        # Save SED figure
        fig_sed.savefig(results_dir + 'sed_plots/' + idx +
                        fsuffix + '_sed.png',
                        bbox_inches='tight')
        plt.close(fig_sed)

        # Save the 500 micron excess dictionary as a pickle
        excess_file = open(results_dir + 'excess500/' + idx +
                           fsuffix + '.pickle', 'wb')
        pickle.dump(excess500, excess_file)
        excess_file.close()

        # Save the MCMC results
        mcmc_results_file = open(results_dir + 'mcmc_fits/' + idx +
                                 fsuffix + '_mcmc.pickle',
                                 'wb')
        pickle.dump(fit_results, mcmc_results_file)
        mcmc_results_file.close()

        print '\tDone!'
        print ''
    print 'Done with all sources!'

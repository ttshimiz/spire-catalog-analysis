# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 20:39:10 2015

Script to run the single temperature blackbody fitting on the BAT AGN that have detections
at all the wavelengths from 160 - 500 microns. However only the 160-350 micron
points will be fit.
After fitting, the 500 micron excess will be measured by comparing the best fit
model 500 micron point with the observed 500 micron point.
The excess will be measured as a percentage in the following way:
    excess500 = [F_obs(500) - F_model(500)] / F_model(500)

The excess values along with fit results will be saved in a text file.

@author: ttshimiz
"""

import sys
sys.path.append('/Users/ttshimiz/Github/bat-agn-sed-fitting/')

import numpy as np
import pandas as pd
import astropy.units as u
import triangle
import matplotlib.pyplot as plt
import seaborn
import pickle
import models as sed_mod
import fitting as sed_fit
import plotting as sed_plot
from astropy.modeling import fitting as apy_fit

seaborn.set()

# Directory where all of the data is stored
data_dir = '/Users/ttshimiz/Github/bat-data/'


# Functions for analysis
def excess(model, observed):
    return (observed - model) / model


def calc_excess(model, obs500):
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
    dummy = model.copy()
    fixed = np.array([dummy.fixed[n] for n in dummy.param_names])
    
    model500 = np.zeros(len(dummy.chain_nb))
    for i in range(len(dummy.chain_nb)):
        dummy.parameters[~fixed] = model.chain_nb[i]
        model500[i] = dummy([500.])

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


def run_all_sources(beta=2.0, nwalkers=50, nburn=200, nsteps=1000,
                    results_dir='results/betaFree/',
                    fsuffix='_betafree_twotemp_500excess'):

    # Upload the Herschel data
    herschel_data = pd.read_csv(data_dir+'bat_herschel.csv', index_col=0)

    # Upload the general information for the BAT sources
    gen_data = pd.read_csv(data_dir+'bat_info.csv', index_col=0)

    # We will use the 160, 250, 350, and 500 micron data
    # We also need the redshift as well
    fluxes = pd.DataFrame({'H70':  herschel_data.PACS70,
                           'H160': herschel_data.PACS160,
                           'H250': herschel_data.PSW,
                           'H350': herschel_data.PMW,
                           'H500': herschel_data.PLW})

    flux_err = pd.DataFrame({'H70_err':  herschel_data.PACS70_err,
                             'H160_err': herschel_data.PACS160_err,
                             'H250_err': herschel_data.PSW_err,
                             'H350_err': herschel_data.PMW_err,
                             'H500_err': herschel_data.PLW_err})

    redshifts = gen_data.Redshift
    lum_dists = gen_data['Dist_[Mpc]']
	
    # Find all of the sources that are detected at every wavelength
    detect_all = (fluxes > 0).all(axis=1)

    # Iterate over all of the sources
    waves = np.array([70., 160., 250., 350., 500.])
    filts = np.array(['PACS70', 'PACS160', 'PSW', 'PMW', 'PLW'])
	
	# Initialize the model and fitter
    base_model = sed_mod.Greybody(0.0, 25., 2.0)
    lev_marq = apy_fit.LevMarLSQFitter()
    bayes = sed_fit.SEDBayesFitter(nwalkers=nwalkers, nsteps=nsteps, nburn=nburn,
                                   threads=4)

	# Fix the emissivity
    base_model.beta.fixed = True
	
#    for idx in fluxes[detect_all].index[2:]:
    for idx in ['ESO244-IG030']:
        print 'Running analysis on ...' + idx
        flux_src = fluxes.loc[idx][['H160', 'H250', 'H350']].values
        flux_err_src = flux_err.loc[idx][['H160_err', 'H250_err', 'H350_err']].values
        zz_src = redshifts.loc[idx]
        ld_src = lum_dists.loc[idx]
        waves_use = np.array([160., 250., 350.])
        filts_use = np.array(['PACS160', 'PSW', 'PMW'])
		
        print '\tFitting model...'
        model_ml = base_model.copy()
        model_ml.set_redshift(zz_src)
        model_ml.set_lumD(ld_src)

        mdust_init = np.log10(flux_src[1]/model_ml([250.]))
        model_ml.mdust.value = mdust_init

        model_init = lev_marq(model_ml, waves_use, flux_src, weights=1/flux_err_src,
                              maxiter=1000)
        model_final = bayes.fit(model_init, waves_use, flux_src, yerr=flux_err_src,
                                filts=filts_use)

        # Calculate the 500 micron excess
        flux500 = fluxes.loc[idx]['H500']
        print '\tCalculating 500 micron excess...'
        excess500 = calc_excess(model_final, flux500)
        
        print '\tPlotting and saving results...'
        
        # Plot a corner triangle figure to look at relationships
        # between parameters
        fig_tri = sed_plot.plot_triangle(model_final)
        
        # Save the triangle plot
        fig_tri.savefig(results_dir + 'triangle_plots/' + idx +
                        fsuffix + '_triangle.png',
                        bbox_inches='tight')
        plt.close(fig_tri)

        # Plot the data along with the median SED. Shade in the 95% confidence
        # interval
        sed_src = fluxes.loc[idx][['H70', 'H160', 'H250', 'H350', 'H500']].values
        sed_src_err = flux_err.loc[idx][['H70_err', 'H160_err', 'H250_err', 'H350_err', 'H500_err']].values
        fig_fit = sed_plot.plot_fit(waves, sed_src, model_final, obs_err=sed_src_err,
                                    plot_components=False, plot_mono_fluxes=True,
                                    filts=filts, plot_fit_spread=True,
                                    name=idx, plot_params=True)
        ax = fig_fit.gca()
        ax.set_xlim(10, 1000)
        
        # Save SED figure
        fig_fit.savefig(results_dir + 'sed_plots/' + idx +
                        fsuffix + '_sed.png',
                        bbox_inches='tight')
        plt.close(fig_fit)

        # Save the 500 micron excess dictionary as a pickle
        excess_file = open(results_dir + 'excess500/' + idx +
                           fsuffix + '.pickle', 'wb')
        pickle.dump(excess500, excess_file)
        excess_file.close()

        # Save the MCMC results
        fit_dict = {'name': idx,
                    'flux': sed_src,
                    'flux_err': sed_src_err,
                    'best_fit_model': model_final,
                    'filters': filts,
                    'waves': waves}

        fit_results_file = open(results_dir + 'fit_results/' + idx +
                                fsuffix + '_fit_result.pickle',
                                'wb')
        pickle.dump(fit_dict, fit_results_file)
        fit_results_file.close()

        print '\tDone!'
        print ''
    print 'Done with all sources!'

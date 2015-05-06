# Script to print out all of the luminosity-luminosity correlation data for ASURV
#There will be six columns in the files, 1st column is independent variable, second column is censor for independent variable,
# third column is dependent variable, fourth is censor indicator for dependent variable,
# fifth is test variable and sixth is censor indicator for test variable
#For the censor indicators:
#    1 : Detected Point
#    0 : Upper limit

from pylab import *

github_dir = '/Users/ttshimiz/Github/spire-catalog-analysis/survival-analysis/'

wavebands = ['70', '160', '250', '350', '500', 'BAT']

for n in range(len(wavebands)):
	for j in range(len(wavebands)):
		if n < j:
			
			execfile('/Users/ttshimiz/Dropbox/Research/Thesis/scripts/upload_bat_ir_database.py')
			if n == 0:
				indep_var = l70
				indep_var_err = l70_err
			elif n == 1:
				indep_var = l160
				indep_var_err = l160_err
			elif n == 2:
				indep_var = l250
				indep_var_err = l250_err
			elif n == 3:
				indep_var = l350
				indep_var_err = l350_err
			elif n == 4:
				indep_var = l500
				indep_var_err = l500_err
			elif n == 5:
				indep_var = lbat
				indep_var_err = ones(len(lbat))
		
			if j == 0:
				dep_var = l70
				dep_var_err = l70_err	
			elif j == 1:
				dep_var = l160
				dep_var_err = l160_err
			elif j == 2:
				dep_var = l250
				dep_var_err = l250_err
			elif j == 3:
				dep_var = l350
				dep_var_err = l350_err
			elif j == 4:
				dep_var = l500
				dep_var_err = l500_err
			elif j == 5:
				dep_var = lbat
				dep_var_err = ones(len(lbat))

			dist_sq = dist_mpc**2
		
			censor_indep = ones(len(indep_var))
			censor_dep = ones(len(indep_var))
			censor_dist = ones(len(indep_var))
		
			undetected_indep = indep_var == 0
			undetected_dep = dep_var == 0
		
			censor_indep[undetected_indep] = 0
			censor_dep[undetected_dep] = 0
		
			indep_var[undetected_indep] = indep_var_err[undetected_indep]
			dep_var[undetected_dep] = dep_var_err[undetected_dep]
			
			#Get rid of radio loud objects
			spire_color_250_350 = h250/h350
			spire_color_350_500 = h350/h500
			a = (spire_color_250_350 < 1.75) & (isfinite(spire_color_250_350)) & (spire_color_250_350 != 0)
			b = (spire_color_350_500 < 1.5) & (isfinite(spire_color_350_500)) & (spire_color_350_500 != 0)
			c = a & b
			
			idx = (indep_var != 0) & (dep_var != 0)
			indep_var = indep_var[idx]
			dep_var = dep_var[idx]
			censor_indep = censor_indep[idx]
			censor_dep = censor_dep[idx]
			censor_dist = censor_dist[idx]
			dist_sq = dist_sq[idx]
		
			asurv_dat = transpose(vstack([log10(indep_var), censor_indep, log10(dep_var), censor_dep, dist_sq, censor_dist]))
			savetxt(github_dir+'/l'+wavebands[n]+'_l'+str(wavebands[j])+'.dat', asurv_dat, fmt = ['%10.6g', '%i', '%10.6g', '%i', '%10.6g', '%i'])

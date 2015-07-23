# Script to print out all of the luminosity-luminosity correlation data for ASURV
#There will be six columns in the files, 1st column is independent variable, second column is censor for independent variable,
# third column is dependent variable, fourth is censor indicator for dependent variable,
# fifth is test variable and sixth is censor indicator for test variable
#For the censor indicators:
#    1 : Detected Point
#    0 : Upper limit

from pylab import *

github_dir = '/Users/ttshimiz/Github/spire-catalog-analysis/survival-analysis/'

dir = '/Users/ttshimiz/Github/bat-data/'
bat_herschel = pd.read_csv(dir+'bat_herschel.csv', index_col=0)
bat_info = pd.read_csv(dir+'bat_info.csv', index_col=0)
bat_flux = pd.read_csv(dir+'bat_bat_flux.csv', index_col=0)

dist_mpc = bat_info['Dist_[Mpc]'].values
h70 = bat_herschel['PACS70'].values
h70_err = bat_herschel['PACS70_err'].values
h160 = bat_herschel['PACS160'].values
h160_err = bat_herschel['PACS160_err'].values
h250 = bat_herschel['PSW'].values
h250_err = bat_herschel['PSW_err'].values
h250_flag = bat_herschel['PSW_flag'].values
h350 = bat_herschel['PMW'].values
h350_err = bat_herschel['PMW_err'].values
h350_flag = bat_herschel['PMW_flag'].values
h500 = bat_herschel['PLW'].values
h500_err = bat_herschel['PLW_err'].values
h500_flag = bat_herschel['PLW_flag'].values

sy1_sample = np.zeros(313, dtype = 'bool')
sy2_sample = np.zeros(313, dtype = 'bool')
liner_sample = np.zeros(313, dtype = 'bool')
agn_sample = np.zeros(313, dtype = 'bool')


#Separate the Sy 1 and Sy 2 samples. Sy 1 = Sy 1, 1.2, and 1.5. Sy 2 = Sy 1.8, 1.9, and 2
#Also create a LINER sample and AGN sample

for i in range(313):

    type_split = bat_info.ix[i, 'BAT_Type'].split()

    if (type_split[0] == 'Sy'):
        if ((type_split[1] == '1') | (type_split[1] == '1.2') | (type_split[1] == '1.4') | (type_split[1] == '1.5')):

            sy1_sample[i] = 1

        elif ((type_split[1] == '1.8') | (type_split[1] == '1.9') | (type_split[1] == '2')):
            sy2_sample[i] = 1
    else:

        if (type_split[0] == 'LINER'):
            liner_sample[i] = 1

        elif (type_split[0] == 'AGN'):
            agn_sample[i] = 1

lbat = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*bat_flux['BAT_flux'].values
l70  = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/70.0e-4)*h70*10**(-23)
l70_err  = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/70.0e-4)*h70_err*10**(-23)
l160 = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/160.0e-4)*h160*10**(-23)
l160_err = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/160.0e-4)*h160_err*10**(-23)          
l250 = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/250.0e-4)*h250*10**(-23)
l250_err = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/250.0e-4)*h250_err*10**(-23)
l350 = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/350.0e-4)*h350*10**(-23)
l350_err = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/350.0e-4)*h350_err*10**(-23)
l500 = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/500.0e-4)*h500*10**(-23)
l500_err = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/500.0e-4)*h500_err*10**(-23)

wavebands = ['70', '160', '250', '350', '500', 'BAT']

for n in range(len(wavebands)):
	for j in range(len(wavebands)):
		if n < j:
			
			#execfile('/Users/ttshimiz/Dropbox/Research/Thesis/scripts/upload_bat_ir_database.py')
			if n == 0:
				indep_var = l70.copy()
				indep_var_err = l70_err.copy()
				indep_var_flag = np.ones(len(indep_var), dtype = 'str')
			elif n == 1:
				indep_var = l160.copy()
				indep_var_err = l160_err.copy()
				indep_var_flag = np.ones(len(indep_var), dtype = 'str')
			elif n == 2:
				indep_var = l250.copy()
				indep_var_err = l250_err.copy()
				indep_var_flag = h250_flag.copy()
			elif n == 3:
				indep_var = l350.copy()
				indep_var_err = l350_err.copy()
				indep_var_flag = h350_flag.copy()
			elif n == 4:
				indep_var = l500.copy()
				indep_var_err = l500_err.copy()
				indep_var_flag = h500_flag.copy()
			elif n == 5:
				indep_var = lbat.copy()
				indep_var_err = ones(len(lbat))
				indep_var_flag = np.ones(len(indep_var), dtype = 'str')
		
			if j == 0:
				dep_var = l70.copy()
				dep_var_err = l70_err.copy()
				dep_var_flag = np.ones(len(indep_var), dtype = 'str')	
			elif j == 1:
				dep_var = l160.copy()
				dep_var_err = l160_err.copy()
				dep_var_flag = np.ones(len(indep_var), dtype = 'str')	
			elif j == 2:
				dep_var = l250.copy()
				dep_var_err = l250_err.copy()
				dep_var_flag = h250_flag.copy()
			elif j == 3:
				dep_var = l350.copy()
				dep_var_err = l350_err.copy()
				dep_var_flag = h350_flag.copy()
			elif j == 4:
				dep_var = l500.copy()
				dep_var_err = l500_err.copy()
				dep_var_flag = h500_flag.copy()
			elif j == 5:
				dep_var = lbat.copy()
				dep_var_err = ones(len(lbat))
				dep_var_flag = np.ones(len(indep_var), dtype = 'str')	

			dist_sq = dist_mpc.copy()**2
		
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
			#spire_color_250_350 = h250/h350
			#spire_color_350_500 = h350/h500
			#a = (spire_color_250_350 < 1.75) & (isfinite(spire_color_250_350)) & (spire_color_250_350 != 0)
			#b = (spire_color_350_500 < 1.5) & (isfinite(spire_color_350_500)) & (spire_color_350_500 != 0)
			#c = a & b
			flag_indep = (indep_var_flag != 'UC') & (indep_var_flag != 'UD') & (indep_var_flag != 'AC') & (indep_var_flag != 'AD')
			flag_dep = (dep_var_flag != 'UC') & (dep_var_flag != 'UD') & (dep_var_flag != 'AC') & (dep_var_flag != 'AD')
			idx = (indep_var != 0) & (dep_var != 0) & flag_indep & flag_dep & sy2_sample
			indep_var = indep_var[idx]
			dep_var = dep_var[idx]
			censor_indep = censor_indep[idx]
			censor_dep = censor_dep[idx]
			censor_dist = censor_dist[idx]
			dist_sq = dist_sq[idx]
		
			asurv_dat = transpose(vstack([log10(indep_var), censor_indep, log10(dep_var), censor_dep, dist_sq, censor_dist]))
			savetxt(github_dir+'/l'+wavebands[n]+'_l'+str(wavebands[j])+'_sy2.dat', asurv_dat, fmt = ['%10.6g', '%i', '%10.6g', '%i', '%10.6g', '%i'])

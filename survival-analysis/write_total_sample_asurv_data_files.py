#Script to write all of the two sample data files to use with ASURV
#The two samples will be Seyfert 1's and 2's
#The format is I1, I2, Data where I1 indicates which group it is in and I2 is the censor indicator
#For the censor indicator: 1 = Lower Limit, 0 = Detection, and -1 = Upper Limit


from pylab import *

github_dir = '/Users/ttshimiz/Github/'

#Upload all of the Herschel data
#execfile('/Users/ttshimiz/Dropbox/Research/Thesis/scripts/upload_bat_ir_database.py')
dir = '/Users/ttshimiz/Github/bat-data/'
bat_herschel = pd.read_csv(dir+'bat_herschel.csv', index_col=0)

h250 = bat_herschel['PSW'].values
h250_err = bat_herschel['PSW_err'].values
h250_flag = bat_herschel['PSW_flag'].values
h350 = bat_herschel['PMW'].values
h350_err = bat_herschel['PMW_err'].values
h350_flag = bat_herschel['PMW_flag'].values
h500 = bat_herschel['PLW'].values
h500_err = bat_herschel['PLW_err'].values
h500_flag = bat_herschel['PLW_flag'].values
      
l250 = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/250.0e-4)*h250*10**(-23)
l250_err = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/250.0e-4)*h250_err*10**(-23)
l350 = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/350.0e-4)*h350*10**(-23)
l350_err = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/350.0e-4)*h350_err*10**(-23)
l500 = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/500.0e-4)*h500*10**(-23)
l500_err = 4*np.pi*(dist_mpc*10**6*3.09e18)**2*(3.0e10/500.0e-4)*h500_err*10**(-23)

wavebands = ['250', '350', '500']

for i in range(len(wavebands)):

	if i == 0:
		cens = np.array(l250 == 0, dtype=int)*-1
		lum = l250
		lum[l250 == 0] = l250_err[l250 == 0]
		flag = (h250_flag != 'UC') & (h250_flag != 'UD') & (h250_flag != 'AC') & (h250_flag != 'AD')
		lum = lum[flag]
		cens = cens[flag]

	elif i == 1:
		cens = np.array(l350 == 0, dtype=int)*-1
		lum = l350
		lum[l350 == 0] = l350_err[l350 == 0]
		flag = (h350_flag != 'UC') & (h350_flag != 'UD') & (h350_flag != 'AC') & (h350_flag != 'AD')		
		lum = lum[flag]
		cens = cens[flag]

	elif i == 2:
		cens = np.array(l500 == 0, dtype=int)*-1
		lum = l500
		lum[l500 == 0] = l500_err[l500 == 0]
		flag = (h500_flag != 'UC') & (h500_flag != 'UD') & (h500_flag != 'AC') & (h500_flag != 'AD')
		lum = lum[flag]
		cens = cens[flag]
			
	one_samp_data  = transpose(vstack([cens, log10(lum)]))
	savetxt(github_dir+'spire-catalog-analysis/survival-analysis/l'+wavebands[i]+'_asurv.dat', one_samp_data, fmt = ['%i', '%10.6g'])
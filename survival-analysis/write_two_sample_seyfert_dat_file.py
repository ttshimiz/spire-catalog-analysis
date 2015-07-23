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

flag_250 = (h250_flag != 'UC') & (h250_flag != 'UD') & (h250_flag != 'AC') & (h250_flag != 'AD')
flag_350 = flag = (h350_flag != 'UC') & (h350_flag != 'UD') & (h350_flag != 'AC') & (h350_flag != 'AD')
flag_500 = flag = (h500_flag != 'UC') & (h500_flag != 'UD') & (h500_flag != 'AC') & (h500_flag != 'AD')

wavebands = ['70', '160', '250', '350', '500']

for i in range(len(wavebands)):

	if i == 0:
		sy1_lum = l70[sy1_sample]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l70_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l70[sy2_sample]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l70_err[sy2_sample][sy2_cens_ind]
		sy1_indicator = zeros(sum(sy1_sample))
		sy2_indicator = ones(sum(sy2_sample))

	elif i == 1:
		sy1_lum = l160[sy1_sample]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l160_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l160[sy2_sample]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l160_err[sy2_sample][sy2_cens_ind]
		sy1_indicator = zeros(sum(sy1_sample))
		sy2_indicator = ones(sum(sy2_sample))

		
	elif i == 2:
		sy1_lum = l250[sy1_sample & flag_250]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l250_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l250[sy2_sample & flag_250]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l250_err[sy2_sample][sy2_cens_ind]
		sy1_indicator = zeros(sum(sy1_sample & flag_250))
		sy2_indicator = ones(sum(sy2_sample & flag_250))

		
	elif i == 3:
		sy1_lum = l350[sy1_sample & flag_350]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l350_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l350[sy2_sample & flag_350]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l350_err[sy2_sample][sy2_cens_ind]
		sy1_indicator = zeros(sum(sy1_sample & flag_350))
		sy2_indicator = ones(sum(sy2_sample & flag_350))

		
	elif i == 4:
		sy1_lum = l500[sy1_sample & flag_500]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l500_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l500[sy2_sample & flag_500]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l500_err[sy2_sample][sy2_cens_ind]
		sy1_indicator = zeros(sum(sy1_sample & flag_500))
		sy2_indicator = ones(sum(sy2_sample & flag_500))
				
	sy1_censor = zeros(len(sy1_indicator))
	sy1_censor[sy1_cens_ind] = -1
	sy2_censor = zeros(len(sy2_indicator))
	sy2_censor[sy2_cens_ind] = -1
	
	two_samp_data  = transpose(vstack([hstack([sy1_indicator, sy2_indicator]), hstack([sy1_censor, sy2_censor]), hstack([log10(sy1_lum), log10(sy2_lum)])]))
	savetxt(github_dir+'spire-catalog-analysis/survival-analysis/two_sample_seyferts_l'+wavebands[i]+'.dat', two_samp_data, fmt = ['%i', '%i', '%10.6g'])

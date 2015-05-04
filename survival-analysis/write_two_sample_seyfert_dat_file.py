#Script to write all of the two sample data files to use with ASURV
#The two samples will be Seyfert 1's and 2's
#The format is I1, I2, Data where I1 indicates which group it is in and I2 is the censor indicator
#For the censor indicator: 1 = Lower Limit, 0 = Detection, and -1 = Upper Limit


from pylab import *

github_dir = '/Users/ttshimiz/Github/'

#Upload all of the Herschel data
execfile('/Users/ttshimiz/Dropbox/Research/Thesis/scripts/upload_bat_ir_database.py')

sy1_indicator = zeros(sum(sy1_sample))
sy2_indicator = ones(sum(sy2_sample))

wavebands = ['70', '160', '250', '350', '500']

for i in range(len(wavebands)):

	if i == 0:
		sy1_lum = l70[sy1_sample]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l70_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l70[sy2_sample]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l70_err[sy2_sample][sy2_cens_ind]
	elif i == 1:
		sy1_lum = l160[sy1_sample]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l160_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l160[sy2_sample]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l160_err[sy2_sample][sy2_cens_ind]
		
	elif i == 2:
		sy1_lum = l250[sy1_sample]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l250_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l250[sy2_sample]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l250_err[sy2_sample][sy2_cens_ind]
		
	elif i == 3:
		sy1_lum = l350[sy1_sample]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l350_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l350[sy2_sample]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l350_err[sy2_sample][sy2_cens_ind]
		
	elif i == 4:
		sy1_lum = l500[sy1_sample]
		sy1_cens_ind = sy1_lum == 0
		sy1_lum[sy1_cens_ind] =  l500_err[sy1_sample][sy1_cens_ind]
		sy2_lum = l500[sy2_sample]
		sy2_cens_ind = sy2_lum == 0
		sy2_lum[sy2_cens_ind] =  l500_err[sy2_sample][sy2_cens_ind]
				
	sy1_censor = zeros(len(sy1_indicator))
	sy1_censor[sy1_cens_ind] = -1
	sy2_censor = zeros(len(sy2_indicator))
	sy2_censor[sy2_cens_ind] = -1
	
	two_samp_data  = transpose(vstack([hstack([sy1_indicator, sy2_indicator]), hstack([sy1_censor, sy2_censor]), hstack([log10(sy1_lum), log10(sy2_lum)])]))
	savetxt(github_dir+'spire-catalog-analysis/survival-analysis/two_sample_seyferts_l'+wavebands[i]+'.dat', two_samp_data, fmt = ['%i', '%i', '%10.6g'])

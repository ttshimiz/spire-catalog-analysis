#Script to write all of the two sample data files to use with ASURV
#The two samples will be Seyfert 1's and 2's
#The format is I1, I2, Data where I1 indicates which group it is in and I2 is the censor indicator
#For the censor indicator: 1 = Lower Limit, 0 = Detection, and -1 = Upper Limit


from pylab import *

github_dir = '/Users/ttshimiz/Github/'

#Upload all of the Herschel data
execfile('/Users/ttshimiz/Dropbox/Research/Thesis/scripts/upload_bat_ir_database.py')

wavebands = ['250', '350', '500']

for i in range(len(wavebands)):

	if i == 0:
		cens = np.array(l250 == 0, dtype=int)*-1
		lum = l250
		lum[l250 == 0] = l250_err[l250 == 0]
	elif i == 1:
		cens = np.array(l350 == 0, dtype=int)*-1
		lum = l350
		lum[l350 == 0] = l350_err[l350 == 0]		
	elif i == 2:
		cens = np.array(l500 == 0, dtype=int)*-1
		lum = l500
		lum[l500 == 0] = l500_err[l500 == 0]
			
	one_samp_data  = transpose(vstack([cens, log10(lum)]))
	savetxt(github_dir+'spire-catalog-analysis/survival-analysis/l'+wavebands[i]+'_asurv.dat', one_samp_data, fmt = ['%i', '%10.6g'])
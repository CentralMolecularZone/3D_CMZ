### Code to calculate various individual cloud extinction method values and store in a new .tex table

#Import libraries
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt 
import matplotlib.cm as cm # colormaps
import matplotlib.colors as colors
import astropy.io.fits as pyfits

###list of file names in order of mask id###

#mising last two clouds because they don't work
cloud_list = ["G359.475-0.044", "G359.508-0.135", "G359.561-0.001", "G359.595-0.223","G359.608+0.018","G359.688-0.132", "G359.701+0.032", "G359.865+0.023", "G359.88-0.081", "G359.979-0.071", "G0.014-0.016","G0.035+0.032", "G0.068-0.076", "G0.105-0.08", "G0.116+0.003","G0.143-0.083", "G0.255+0.02", "G0.327-0.195", "G0.342+0.06",  "G0.342-0.085", "G0.379+0.05", "G0.413+0.048", "G0.488+0.008", "G0.645+0.03",  "G0.666-0.028", "G0.716-0.09", "G0.816-0.185" , "G0.888-0.044", "G1.075-0.049", "G1.601+0.012", "G1.652-0.052"  ]

flux_diff = []
flux_ratio = []
ext_fraction = []
I_model = []
I_cloud = []
min_flux = []
smoothed_flux_ratio = []

###Import Files###
for i in cloud_list:
	Cloud = i


	##import FITS files##

	#Spitz resid
	spitz_fits = pyfits.open('./Cloud_masks/{}/{}_spitzer_cutout_isolated.fits'.format(Cloud, Cloud))

	#Spitz RAW
	spitz_RAW_fits = pyfits.open('./Cloud_masks/{}/{}_spitzer_RAW_cutout_isolated.fits'.format(Cloud, Cloud))

	#Spitz Smoothed
	spitz_smoothed_fits = pyfits.open('./Cloud_masks/{}/{}_smoothed_cutout_isolated.fits'.format(Cloud, Cloud))

	#Herschel
	hers_fits = pyfits.open('./Cloud_masks/{}/{}_herschel_cutout_isolated_regrid_to_spitzer.fits'.format(Cloud, Cloud))





	#Model from 7 to 9kpc, cutout for cloud region
	#model7to9_fits = pyfits.open('./Cloud_masks/{}/{}_model_7to9kpc_isolated_cutout.fits'.format(Cloud, Cloud))


	# copy the FITS data into a numpy array
	#spitz_data = spitz[0].data
	hers = hers_fits[0].data
	spitz = spitz_fits[0].data
	spitz_RAW = spitz_RAW_fits[0].data
	spitz_smoothed = spitz_smoothed_fits[0].data

	#model7to9 = model7to9_fits[0].data


	h_spitz = pyfits.getheader('./Cloud_masks/{}/{}_spitzer_cutout_isolated.fits'.format(Cloud, Cloud))








	###Calculate flux fraction###
	def flux_frac(cloud_spitz, cloud_model,h_cloud_spitz):

		cloud_flux_frac_array = cloud_spitz/cloud_model

		#fraction of pixels that are optically thick (i.e. less spitzer flux compared to the model)
		thick_pixel_frac = np.size(np.where(cloud_flux_frac_array<1))/np.size(np.where(~np.isnan(cloud_spitz) == True))

    
		return thick_pixel_frac, cloud_flux_frac_array

	def flux_ratio_func(cloud_spitz, cloud_model,h_cloud_spitz):

		cloud_flux_ratio_array = (cloud_spitz - 52.09)/cloud_model

		return cloud_flux_ratio_array

	def flux_diff_func(cloud_spitz, cloud_model,h_cloud_spitz):

		cloud_flux_diff_array = (cloud_model - cloud_spitz) - 72.77

		return cloud_flux_diff_array

	extinction_frac, flux_frac_array = flux_frac(spitz, spitz_smoothed, h_spitz)
	flux_ratio_array = flux_ratio_func(spitz,spitz_smoothed, h_spitz)
	flux_diff_array = flux_diff_func(spitz,spitz_smoothed, h_spitz)




	#make FITS file of flux fractions (spitzer/model)
	pyfits.writeto('./Cloud_masks/{}/{}_flux_frac.fits'.format(Cloud, Cloud), flux_frac_array, h_spitz, overwrite=True)


	pyfits.writeto('./Cloud_masks/{}/{}_flux_ratio.fits'.format(Cloud, Cloud), flux_ratio_array, h_spitz, overwrite=True)

	pyfits.writeto('./Cloud_masks/{}/{}_flux_diff.fits'.format(Cloud, Cloud), flux_diff_array, h_spitz, overwrite=True)

	#print("{} Median flux difference = ".format(Cloud), (np.nanmedian(spitz_smoothed)-np.nanmedian(spitz)))

	#print("{} Median Flux Ratio fore corrected = ".format(Cloud), (np.nanmedian(spitz) -49.792)/np.nanmedian(spitz_smoothed))
	

	#print("{} Flux Diff = ".format(Cloud), (np.nanmedian(spitz_smoothed)-np.nanmedian(spitz)) - 72.77)

	#print("{} Min Flux = ".format(Cloud), np.nanmin(spitz_RAW))

	#print("{} I_cloud = ".format(Cloud), np.nanmedian(spitz))

	print("{} Smoothed Flux Ratio =".format(Cloud),(np.nanmedian(spitz_smoothed) - 52.09)/np.nanmedian(spitz_smoothed))
	
	flux_diff.append(( np.nanmedian(spitz_smoothed)-np.nanmedian(spitz)) - 72.77) #subtract (back+fore based on smoothed image)

	# Foreground = min(SgB2) + min(Brick)/2 = 38.64+65.54/2 = 52.09
	
	flux_ratio.append((np.nanmedian(spitz) - 52.09)/np.nanmedian(spitz_smoothed)) #subtract foreground estimate
	
	ext_fraction.append(extinction_frac)
	
	min_flux.append(np.nanmin(spitz_RAW))

	I_cloud.append(np.nanmedian(spitz))
	I_model.append(np.nanmedian(spitz_smoothed))
	
	smoothed_flux_ratio.append((np.nanmedian(spitz_smoothed) - 52.09)/np.nanmedian(spitz_smoothed))



print("Median Smoothed Flux Ratio =", np.nanmedian(smoothed_flux_ratio))
print("Mean Smoothed Flux Ratio =", np.nanmean(smoothed_flux_ratio))

#### Write to .txt or .tex file so I don't have to 
#### keep doing this by hand...

import pandas as pd

CMZ_cat_df = pd.read_csv("CMZ_cloud_catalogue_data.csv", usecols = ["leaf_id","l","b"])




d = {'leaf id':CMZ_cat_df.leaf_id, 'l': CMZ_cat_df.l, 'b':CMZ_cat_df.b,   'I_cloud': I_cloud, 'I_model': I_model,'Flux Ratio': flux_ratio, 'Extinction Fraction': ext_fraction, 'Flux Difference':flux_diff,  'Min Flux':min_flux}

df = pd.DataFrame(data=d)


df.to_csv('CMZ_cat_flux_method_vals_updated.txt', sep='\t', index = False)

df.to_latex('CMZ_cat_flux_method_vals_updated.tex',index=False, float_format='%.3f')





"""
fig, ax = plt.subplots()

hist, hist_edges = np.histogram(flux_diff)
plt.hist(flux_diff)

plt.xlabel("flux difference bins",fontsize=18)
plt.ylabel("# of clouds",fontsize=18)

plt.show()
"""
"""
####Calculate Degree of Extinction Sigma###

#Find sigma value of the flux_frac image

flux_frac_stdv = np.nanstd(flux_frac_array)
#print("flux frac stdv = ", flux_frac_stdv)

#compare median value of flux frac to stdv?
#print("median flux frac = ", np.nanmedian(flux_frac_array))

fig, ax = plt.subplots()

hist, hist_edges = np.histogram(flux_frac_array[~np.isnan(flux_frac_array)])
plt.hist(flux_frac_array[~np.isnan(flux_frac_array)], range = (0,1))
plt.vlines(flux_frac_stdv, ymin= 0, ymax = np.max(hist), color = 'r', linestyle = '--')
#plt.vlines(-flux_frac_stdv, ymin= 0, ymax = np.max(hist), color = 'r', linestyle = '--')
plt.vlines(2*flux_frac_stdv, ymin= 0, ymax = np.max(hist), color = 'g', linestyle = '--')
#plt.vlines(-2*flux_frac_stdv, ymin= 0, ymax = np.max(hist), color = 'g', linestyle = '--')
plt.vlines(3*flux_frac_stdv, ymin= 0, ymax = np.max(hist), color = 'y', linestyle = '--')
#plt.vlines(-3*flux_frac_stdv, ymin= 0, ymax = np.max(hist), color = 'y', linestyle = '--')

plt.vlines(np.nanmedian(flux_frac_array), ymin = 0, ymax=np.max(hist), color = 'black', linestyle ='-')

plt.xlabel("flux ratio bins",fontsize=18)
plt.ylabel("counts",fontsize=18)

plt.show()
"""




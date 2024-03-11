#!/usr/bin/env python
# coding: utf-8



#=============================#
#   Dani Lipman (11.12.2020)
#=============================#




# import needed extensions
import numpy as np
import matplotlib.pyplot as plt # plotting package
import matplotlib.cm as cm # colormaps
import matplotlib.colors as colors
import astropy.io.fits as pyfits

cloud_list=["G0.488+0.008", "G0.255+0.02", "G0.342-0.085", "G359.701+0.032", "G0.888-0.044", "G0.342+0.06", 
		"G0.014-0.016", "G0.035+0.032", "G0,105-0.080", "G0.666-0.028"]


for Cloud in cloud_list:
	##OPEN MAPS##
	#fpeak_sub_file = './Cloud_masks/{}/{}_herschel70um_sources_sub_fpeak_cutout_isolated.fits'.format(Cloud, Cloud)
	#fpeak_sub_data = pyfits.open(fpeak_sub_file)[0].data

	nans_sub_file = './Cloud_masks/{}/{}_herschel70um_sources_to_nans_cutout_isolated.fits'.format(Cloud, Cloud)
	nans_sub_data = pyfits.open(nans_sub_file)[0].data

	hers70um = './new_mask/destripe_l000_blue_wgls_rcal_cropped.fits'





	##Import libraries for smoothing##
	from astropy.convolution import Gaussian2DKernel
	from scipy.signal import convolve as scipy_convolve
	from astropy.convolution import convolve_fft
	from scipy import ndimage, misc


	# Make your Gaussian
	# sigma_gauss = sqrt( resolution you want^2 - resolution you have^2) = gaussian sigma
	# gaussian sigma / number of arcseconds per pixel = 2 pixels # Your number will be different


	#Want to match 36" convolution, hers 70um resolution is 6"
	sig_gauss = np.sqrt(36**2 - 6**2 )
	kern = sig_gauss/(3.2) #arcsec/pix

	#h_fpeak = pyfits.getheader(fpeak_sub_file)
	h_nans = pyfits.getheader(nans_sub_file)
	h_hers70 = pyfits.getheader(hers70um)



	#smoothed_fpeak = convolve_fft(fpeak_sub_data, Gaussian2DKernel(x_stddev=kern), normalize_kernel=True, preserve_nan=False, allow_huge=True)

	smoothed_nans = convolve_fft(nans_sub_data, Gaussian2DKernel(x_stddev=kern), normalize_kernel=True, preserve_nan=False, allow_huge=True)

	#pyfits.writeto('./Cloud_masks/{}/{}_herschel70um_sources_sub_fpeak_cutout_smoothed_conv36.fits'.format(Cloud,Cloud), smoothed_fpeak, h_fpeak, overwrite=True) 

	pyfits.writeto('./Cloud_masks/{}/{}_herschel70um_sources_to_nans_cutout_smoothed_conv36.fits'.format(Cloud,Cloud), smoothed_nans, h_nans, overwrite=True) 









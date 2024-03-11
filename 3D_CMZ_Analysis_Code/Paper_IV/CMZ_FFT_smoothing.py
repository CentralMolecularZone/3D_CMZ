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



##OPEN MAPS##

spitz = './COMBINED_GLM/Resid_Mosaic_crop_smooth/GLM_resid_I4_mosaic_cutout_masked_full.fits'
hers70um = './new_mask/destripe_l000_blue_wgls_rcal_cropped_masked_full.fits'
hers70um_l2 = './HIGAL70_l2/destripe_l002_blue_wgls_rcal_cropped_small_masked_full.fits'





##Import libraries for smoothing##
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve_fft
from scipy import ndimage, misc


###SPITZER BACKGROUND SMOOTHING####


background_file = pyfits.open(spitz) 

## copy the FITS data into a numpy array##
background = background_file[0].data
h_spitz = pyfits.getheader(spitz)

# Make your Gaussian
# sigma_gauss = sqrt( resolution you want^2 - resolution you have^2) = gaussian sigma
# gaussian sigma / number of arcseconds per pixel = 2 pixels # Your number will be different
r_eff = 3 #arcmin
res = r_eff/np.sqrt(2)
test_res = res * 60 #arcmin * (60"/arcmin), want in arcseconds

# sigma = FWHM/2.355
spitz_res = 1.98/2.355 #8um channel FWHM=1.98 

sig_gauss = np.sqrt(test_res**2 - spitz_res**2 )
kern = sig_gauss/(1.2) #arcsec/pix



smoothed = convolve_fft(background, Gaussian2DKernel(x_stddev=kern), normalize_kernel=True, preserve_nan=False, allow_huge=True)


pyfits.writeto('./SMOOTHED_MAPS/GLM_resid_I4_mosaic_cutout_FFT_smoothed_3arcmin_normed_kernel.fits', smoothed, h_spitz, overwrite=True) 




####HERSCHEL 70um L0 SMOOTHING####

background_file = pyfits.open(hers70um) 

## copy the FITS data into a numpy array##
background = background_file[0].data
h_hers70 = pyfits.getheader(hers70um)

r_eff = 3 #arcmin
res = r_eff/np.sqrt(2)
test_res = res * 60 #arcmin * (60"/arcmin), want in arcseconds

# sigma = FWHM/2.355
hers70_res = 6/2.355 #70um channel FWHM=6"

sig_gauss = np.sqrt(test_res**2 - hers70_res**2 )
kern = sig_gauss/(3.2) #arcsec/pix



smoothed = convolve_fft(background, Gaussian2DKernel(x_stddev=kern), normalize_kernel=True, preserve_nan=False, allow_huge=True)


pyfits.writeto('./SMOOTHED_MAPS/herschel70um_cutout_FFT_smoothed_3arcmin_normed_kernel.fits', smoothed, h_hers70, overwrite=True) 











####HERSCHEL 70um L2 SMOOTHING####

background_file = pyfits.open(hers70um_l2) 

## copy the FITS data into a numpy array##
background = background_file[0].data
h_hers70_l2 = pyfits.getheader(hers70um_l2)

r_eff = 3 #arcmin
res = r_eff/np.sqrt(2)
test_res = res * 60 #arcmin * (60"/arcmin), want in arcseconds

# sigma = FWHM/2.355
hers70_res = 6/2.355 #70um channel FWHM=6" 


sig_gauss = np.sqrt(test_res**2 - hers70_res**2 )
kern = sig_gauss/(3.2) #arcsec/pix



smoothed = convolve_fft(background, Gaussian2DKernel(x_stddev=kern), normalize_kernel=True, preserve_nan=False, allow_huge=True)


pyfits.writeto('./SMOOTHED_MAPS/destripe_l002_blue_wgls_rcal_cropped_small_masked_full_SMOOTHED.fits', smoothed, h_hers70_l2, overwrite=True) 





#plt.figure(2)
#plt.clf()

#fig2, ax2 = plt.subplots(1,1,dpi=300)
#ax2.imshow(smoothed, cmap=colmap, aspect='auto')

#plt.savefig('smoothed.png'.format(res,contour))





#plotting do-dads 
"""
plt.ion() # do plots in interactive mode    
colmap = plt.get_cmap('plasma') # load gray colormap

# plot the image on the screen
plt.figure(1)
plt.clf()

fig1, ax1 = plt.subplots(1,1,figsize=(3,5),dpi=300)
ax1.imshow(img, cmap=colmap, aspect='auto')

"""






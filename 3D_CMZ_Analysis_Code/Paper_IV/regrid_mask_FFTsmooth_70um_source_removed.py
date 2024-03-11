### This is a python script for use with the Common Astronomical Software Application (CASA)
### that was used to regrid the individual cloud masks and apply them to the Herschel and Spitzer maps

# import needed extensions
import os
import numpy as np
import pyfits


cloud_list = ["G359.475-0.044", "G359.508-0.135", "G359.561-0.001", "G359.595-0.223","G359.608+0.018","G359.688-0.132", "G359.701+0.032", "G359.865+0.023", "G359.88-0.081", "G359.979-0.071", "G0.014-0.016","G0.035+0.032", "G0.068-0.076", "G0.105-0.08", "G0.116+0.003","G0.143-0.083", "G0.255+0.02", "G0.327-0.195", "G0.342+0.06",  "G0.342-0.085", "G0.379+0.05", "G0.413+0.048", "G0.488+0.008", "G0.645+0.03",  "G0.666-0.028", "G0.716-0.09", "G0.816-0.185" , "G0.888-0.044", "G1.075-0.049", "G1.601+0.012", "G1.652-0.052"  ]

cloud_list=["G0.255+0.02"]

ffore_list = np.linspace(0.4,0.6,14)
#ffore_list=[0.50]

#cloud_list=["G1.601+0.012", "G1.652-0.052" ]


###FITS names###


for cloud in cloud_list:
	os.chdir('/mnt/data0/home/lipman/3D_CMZ/Cloud_masks')
	os.chdir(cloud)

	if len(cloud)>12:	
		mask_name = "{}_{}".format(cloud[:7],cloud[9:])
	if len(cloud)<13:
		mask_name = "{}_{}".format(cloud[:6],cloud[8:])


	herschel_image = "column_properunits_conv36_source_only"
	smoothed = "../../SMOOTHED_MAPS/GLM_resid_I4_mosaic_cutout_FFT_smoothed_3arcmin_normed_kernel"

	herschel70um = "../../new_mask/destripe_l000_blue_wgls_rcal"
	herschel70um_smoothed='../../SMOOTHED_MAPS/herschel70um_cutout_FFT_smoothed_3arcmin_normed_kernel'



	####Now for NANs Image####
	

	for ffore in ffore_list:
		####Import the Convolved to_nans map###
		importfits(fitsimage='{}_ExtN70um_sources_to_nans_ffore{:.2f}_cutout_smoothed_conv36.fits'.format(cloud,ffore), 
			imagename='{}_ExtN70um_sources_to_nans_ffore{:.2f}_cutout_smoothed_conv36.im'.format(cloud, ffore),
			overwrite=True
			)

		
		####Regrid Convolved map to Herschel HiGAL###
		imregrid(imagename='{}_ExtN70um_sources_to_nans_ffore{:.2f}_cutout_smoothed_conv36.im'.format(cloud,ffore),
			template = '../../new_mask/column_properunits_conv36_source_only.im',
			output='{}_ExtN70um_sources_to_nans_ffore{:.2f}_cutout_smoothed_conv36_regrid.im'.format(cloud,ffore),
			overwrite=True
		)

		###imsubimage regridded conv map using herschel mask###
		imsubimage(imagename='{}_ExtN70um_sources_to_nans_ffore{:.2f}_cutout_smoothed_conv36_regrid.im'.format(cloud,ffore),
			outfile='{}_ExtN70um_sources_to_nans_ffore{:.2f}_cutout_smoothed_conv36_regrid_isolated.im'.format(cloud,ffore),
			region='{}_cutout_region.crtf'.format(cloud),
			overwrite=True)

		###Export image###
		exportfits(imagename='{}_ExtN70um_sources_to_nans_ffore{:.2f}_cutout_smoothed_conv36_regrid_isolated.im'.format(cloud,ffore),
			fitsimage='{}_ExtN70um_sources_to_nans_ffore{:.2f}_cutout_smoothed_conv36_regrid_isolated.fits'.format(cloud,ffore),
			overwrite=True)












	



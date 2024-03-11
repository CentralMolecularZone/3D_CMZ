### This is a python script for use with the Common Astronomical Software Application (CASA)
### that was used to regrid the individual cloud masks and apply them to the Herschel and Spitzer maps

###Python libraries###
import numpy as np
import os

cloud_list = ["G359.475-0.044", "G359.508-0.135", "G359.561-0.001", "G359.595-0.223","G359.608+0.018","G359.688-0.132", "G359.701+0.032", "G359.865+0.023", "G359.88-0.081", "G359.979-0.071", "G0.014-0.016","G0.035+0.032", "G0.068-0.076", "G0.105-0.08", "G0.116+0.003","G0.143-0.083", "G0.255+0.02", "G0.327-0.195", "G0.342+0.06",  "G0.342-0.085", "G0.379+0.05", "G0.413+0.048", "G0.488+0.008", "G0.645+0.03",  "G0.666-0.028", "G0.716-0.09", "G0.816-0.185" , "G0.888-0.044", "G1.075-0.049" ]



#cloud_list = ["G1.601+0.012", "G1.652-0.052"]

#cloud_list=["G1.075-0.049"]

###FITS names###


for cloud in cloud_list:
	os.chdir('/mnt/data0/home/lipman/3D_CMZ/Cloud_masks')
	os.chdir(cloud)

	if len(cloud)>12:	
		mask_name = "{}_{}".format(cloud[:7],cloud[9:])
	if len(cloud)<13:
		mask_name = "{}_{}".format(cloud[:6],cloud[8:])


	spitz_raw = "../../COMBINED_GLM/Raw_Mosaic_crop/GLM_raw_I4_mosaic_cropped"

	#model_FULL = "../../gaussian_model/gaussian_model_dist_fudgedtomodel_xstdv_lbr_655x590_7to9kpc_regrid"

	#spitzer_image = "../../new_mask/0001_267.33439000_-27.56022000_GLM_00300+0000_mosaic_I4"

	spitzer_image = "../../COMBINED_GLM/Resid_Mosaic_crop_smooth/GLM_resid_I4_mosaic_cutout"

	herschel_image = "column_properunits_conv36_source_only"

	smoothed = "../../SMOOTHED_MAPS/GLM_resid_I4_mosaic_cutout_FFT_smoothed_3arcmin_normed_kernel"

	herschel70um = "../../new_mask/destripe_l000_blue_wgls_rcal"
	#herschel70um = "../../HIGAL70_l2/destripe_l002_blue_wgls_rcal_cropped_small_masked_full"

	herschel70um_smoothed='../../SMOOTHED_MAPS/herschel70um_cutout_FFT_smoothed_3arcmin_normed_kernel'
	#herschel70um_smoothed="../../SMOOTHED_MAPS/destripe_l002_blue_wgls_rcal_cropped_small_masked_full_SMOOTHED"

	herschel70um_sources_to_nans ="../../70um_sources_to_nans"
	#herschel70um_sources_to_nans ="../../70um_sources_to_nans_l2"



	#near_map = "gaussian_model_dust_ridge_e_f_foreground_7to7.911kpc"
	#far_map = "gaussian_model_dust_ridge_e_f_foreground_7to8.089kpc"


	#Import Raw Spitz map
	print("Import Raw Spitz map")
	importfits(fitsimage="{}.fits".format(spitz_raw),
			imagename="{}.im".format(spitz_raw), 
			overwrite = True
				)  

	print(" ")
	print(" ")

	###imsubimage to cut out mask###
	print("imsubimage to cut out mask")
	importfits(fitsimage="{}_mask.fits".format(cloud),
			imagename="{}_mask.im".format(cloud), 
			overwrite = True
				)  
	#CASA's "mask" keyword cannot understand -/+ in file names!!!
	imsubimage(imagename="{}_mask.im".format(cloud), 	  
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_cutout.im".format(mask_name), 
		overwrite = True
			)   

	exportfits(imagename="{}_cutout.im".format(mask_name),
			fitsimage="{}_cutout.fits".format(mask_name), 
			overwrite = True
				)  

	print(" ")
	print(" ")

	###Import mask and regrid to smoothed and raw Spitzer map###
	print("Import mask and regrid to smoothed and raw Spitzer map")

	print(" ")
	print(" ")

	#SMOOTHED
	print("SMOOTHED: Import mask and regrid to smoothed and raw Spitzer map")
	importfits(fitsimage="{}_cutout.fits".format(mask_name),
			imagename="{}_cutout.im".format(mask_name), 
			overwrite = True
				)  

	imregrid(imagename="{}_cutout.im".format(mask_name), 
				template="{}.im".format(smoothed),
				output="{}_cutout_regrid.im".format(mask_name),
				overwrite = True
					)

	exportfits(imagename="{}_cutout_regrid.im".format(mask_name),
			fitsimage="{}_cutout_regrid.fits".format(mask_name),
			overwrite=True
				)
	print(" ")
	print(" ")


	#RAW 
	print("RAW: Import mask and regrid to smoothed and raw Spitzer map")
	imregrid(imagename="{}_cutout.im".format(mask_name), 
				template="{}.im".format(spitz_raw),
				output="{}_cutout_regrid_RAW.im".format(mask_name),
				overwrite = True
					)

	exportfits(imagename="{}_cutout_regrid.im".format(mask_name),
			fitsimage="{}_cutout_regrid_RAW.fits".format(mask_name),
			overwrite=True
				)

	print(" ")
	print(" ")


	###Manually create region file for cloud mask###
	print("Did you manually create the .crtf casa region file??")

	print(" ")
	print(" ")

	###imsubimage smoothed Spitzer image with regridded cloud mask###
	print("imsubimage smoothed Spitzer image with regridded cloud mask")
	imsubimage(imagename="{}.im".format(smoothed), 	  
		mask="{}_cutout_regrid.im".format(mask_name), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_smoothed_cutout_isolated.im".format(cloud),
		overwrite = True
			)                  

	exportfits(imagename="{}_smoothed_cutout_isolated.im".format(cloud),
			fitsimage="{}_smoothed_cutout_isolated.fits".format(cloud),
			overwrite=True
				)


	print(" ")
	print(" ")

	###imsubimage Spitzer resid and raw images with mask###


	#RESID
	print("RESID: imsubimage Spitzer resid and raw images with mask")
	importfits(fitsimage="{}.fits".format(spitzer_image),
			imagename="{}.im".format(spitzer_image), 
			overwrite = True
				)  

	imsubimage(imagename="{}.im".format(spitzer_image), 	  
		mask="{}_cutout_regrid.im".format(mask_name), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_spitzer_cutout_isolated.im".format(cloud),
		overwrite = True
			)                  

	exportfits(imagename="{}_spitzer_cutout_isolated.im".format(cloud),
			fitsimage="{}_spitzer_cutout_isolated.fits".format(cloud),
			overwrite=True
				)
	print(" ")
	print(" ")

	#RAW 
	print("RAW: imsubimage Spitzer resid and raw images with mask")
	imsubimage(imagename="{}.im".format(spitz_raw), 	  
		mask="{}_cutout_regrid_RAW.im".format(mask_name), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_spitzer_RAW_cutout_isolated.im".format(cloud),
		overwrite = True
			)                  

	exportfits(imagename="{}_spitzer_RAW_cutout_isolated.im".format(cloud),
			fitsimage="{}_spitzer_RAW_cutout_isolated.fits".format(cloud),
			overwrite=True
				)

	print(" ")
	print(" ")



	###imregrid cloud mask to herschel; will then apply imsubimage with herschel-gridded cloud mask###
	print("imregrid cloud mask to herschel col den")
	importfits(fitsimage="../../new_mask/{}.fits".format(herschel_image),
			imagename="../../new_mask/{}.im".format(herschel_image), 
			overwrite = True
				)  

	imregrid(imagename="{}_cutout.im".format(mask_name), 
				template="../../new_mask/{}.im".format(herschel_image),
				output="{}_cutout_regrid_to_herschel.im".format(mask_name),
				overwrite = True
					)

	exportfits(imagename = "{}_cutout_regrid_to_herschel.im".format(mask_name),
				fitsimage = "{}_cutout_regrid_to_herschel.fits".format(mask_name),
				overwrite = True
					)

	print(" ")
	print(" ")

	###imregrid cloud mask to herschel 70um; will then apply imsubimage with herschel-gridded cloud mask###
	print("imregrid cloud mask to herschel 70um")

	imregrid(imagename="{}_cutout.im".format(mask_name), 
				template="{}.im".format(herschel70um),
				output="{}_cutout_regrid_to_herschel70um.im".format(mask_name),
				overwrite = True
					)

	exportfits(imagename = "{}_cutout_regrid_to_herschel70um.im".format(mask_name),
				fitsimage = "{}_cutout_regrid_to_herschel70um.fits".format(mask_name),
				overwrite = True
					)

	print(" ")
	print(" ")

	###imsubimage herschel image with herschel-regridded mask###
	print("imsubimage herschel image with herschel-regridded mask")
	imsubimage(imagename="../../new_mask/{}.im".format(herschel_image), 	  
		mask="{}_cutout_regrid_to_herschel.im".format(mask_name), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_{}_cutout.im".format(cloud,herschel_image),
		overwrite = True
			)                  

	exportfits(imagename="{}_{}_cutout.im".format(cloud,herschel_image),
			fitsimage="{}_herschel_cutout_isolated.fits".format(cloud),
			overwrite=True
				) #this file has the herschel pixel scaling and must be regridded to spitzer


	print(" ")
	print(" ")


	#HERSCHEL 70um 
	print("70um: imsubimage Herschel 70um with regridded mask")
	imsubimage(imagename="{}.im".format(herschel70um), 	  
		mask="{}_cutout_regrid_to_herschel70um.im".format(mask_name), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_herschel70um_cutout_isolated.im".format(cloud),
		overwrite = True
			)                  

	exportfits(imagename="{}_herschel70um_cutout_isolated.im".format(cloud),
			fitsimage="{}_herschel70um_cutout_isolated.fits".format(cloud),
			overwrite=True
				)

	print(" ")
	print(" ")

	#70um SMOOTHED
	print("70um: imregrid cloud mask to Herschel 70um")
	imregrid(imagename="{}_cutout.im".format(mask_name), 
				template="{}.im".format(herschel70um_smoothed),
				output="{}_cutout_regrid_to_herschel70um_smoothed.im".format(mask_name),
				overwrite = True
					)

	print(" ")
	print(" ")

	print("70um: imsubimage SMOOTHED Herschel 70um with regridded mask")
	imsubimage(imagename="{}.im".format(herschel70um_smoothed), 	  
		mask="{}_cutout_regrid_to_herschel70um_smoothed.im".format(mask_name), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_herschel70um_SMOOTHED_cutout_isolated.im".format(cloud),
		overwrite = True
			)                  

	exportfits(imagename="{}_herschel70um_SMOOTHED_cutout_isolated.im".format(cloud),
			fitsimage="{}_herschel70um_SMOOTHED_cutout_isolated.fits".format(cloud),
			overwrite=True
				)

	print(" ")
	print(" ")

	###imregrid masked-herschel cutout to spitzer###
	print("imregrid masked-herschel cutout to spitzer")
	imregrid(imagename="{}_{}_cutout.im".format(cloud,herschel_image), 
				template="{}.im".format(spitzer_image),
				output="{}_herschel_cutout_regrid_to_spitzer.im".format(cloud),
				overwrite = True
					) #this file is regridded, but zoomed out. Use imsubimage for a region cutout

	imsubimage(imagename = "{}_herschel_cutout_regrid_to_spitzer.im".format(cloud),
				region = "{}_cutout_region.crtf".format(cloud),
				outfile = "{}_herschel_cutout_isolated_regrid_to_spitzer.im".format(cloud),
				overwrite = True
					) #cuts out the smaller region of the regridded herschel masked cloud

	exportfits(imagename="{}_herschel_cutout_isolated_regrid_to_spitzer.im".format(cloud),
			fitsimage="{}_herschel_cutout_isolated_regrid_to_spitzer.fits".format(cloud),
			overwrite=True
				)

	print(" ")
	print(" ")

	###imregrid masked-herschel cutout to 70um###
	print("imregrid masked-herschel cutout to 70um")
	imregrid(imagename="{}_{}_cutout.im".format(cloud,herschel_image), 
				template="{}.im".format(herschel70um),
				output="{}_herschel_cutout_regrid_to_70um.im".format(cloud),
				overwrite = True
					) #this file is regridded, but zoomed out. Use imsubimage for a region cutout

	imsubimage(imagename = "{}_herschel_cutout_regrid_to_70um.im".format(cloud),
				region = "{}_cutout_region.crtf".format(cloud),
				outfile = "{}_herschel_cutout_isolated_regrid_to_70um.im".format(cloud),
				overwrite = True
					) #cuts out the smaller region of the regridded herschel masked cloud

	exportfits(imagename="{}_herschel_cutout_isolated_regrid_to_70um.im".format(cloud),
			fitsimage="{}_herschel_cutout_isolated_regrid_to_70um.fits".format(cloud),
			overwrite=True
				)

	print(" ")
	print(" ")










	####LOAD IN HERSCHEL SUBTRACTED MAPS####


	#TO NANS
	importfits(fitsimage="{}.fits".format(herschel70um_sources_to_nans),
			imagename = "{}.im".format(herschel70um_sources_to_nans),
			overwrite = True)


	####REGRID MASK TO CROPPED HERSCHEL 70UM####
	imregrid(imagename="{}_cutout.im".format(mask_name), 
				template="{}.im".format(herschel70um_sources_to_nans),
				output="{}_cutout_regrid_to_herschel70um_cropped.im".format(mask_name),
				overwrite = True
					)


	####CUTOUT CLOUDS FROM SUB/NAN MAPS USING 70UM REGRIDDED MASK####

             

	#TO NANS
	print("70um: imsubimage Herschel 70um TO NANS with regridded mask")
	imsubimage(imagename="{}.im".format(herschel70um_sources_to_nans), 	  
		mask="{}_cutout_regrid_to_herschel70um_cropped.im".format(mask_name), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_herschel70um_sources_to_nans_cutout_isolated.im".format(cloud),
		overwrite = True
			)    



	#70um NANS subtracted maps
	exportfits(imagename="{}_herschel70um_sources_to_nans_cutout_isolated.im".format(cloud),
			fitsimage="{}_herschel70um_sources_to_nans_cutout_isolated.fits".format(cloud),
			overwrite=True
				)

	print(" ")
	print(" ")


	"""

	###imsubimage model with mask###

	###First, we need to integrate the gaussian model along the line of sight to 
	#an approximate distance of the cloud --> use savefits.py###



	###import near/far integrated model maps###

	#import NEAR map
	importfits(fitsimage="./{}.fits".format(near_map),
			imagename="./{}.im".format(near_map), 
			overwrite = True
				)  

	#import FAR map
	importfits(fitsimage="./{}.fits".format(far_map),
			imagename="./{}.im".format(far_map), 
			overwrite = True
				)  

	###regrid near/far model maps to smoothed image###

	#NEAR map regrid to smoothed
	imregrid(imagename="./{}.im".format(near_map), 
				template="../../{}.im".format(smoothed),
				output="./{}_regrid.im".format(near_map),
				overwrite = True
					)
	#FAR map regrid to smoothed
	imregrid(imagename="./{}.im".format(far_map), 
				template="../../{}.im".format(smoothed),
				output="./{}_regrid.im".format(far_map),
				overwrite = True
					)

	#regridded near/far maps export to fits
	exportfits(imagename="./{}_regrid.im".format(near_map),
			fitsimage="./{}_regrid.fits".format(near_map),
			overwrite=True
				)

	exportfits(imagename="./{}_regrid.im".format(far_map),
			fitsimage="./{}_regrid.fits".format(far_map),
			overwrite=True
				)



	###imsubimage near/far model maps to cloud mask and region###

	#imsubimage NEAR map
	imsubimage(imagename="./{}_regrid.im".format(near_map), 	  
		mask="{}_cutout_regrid.im".format(cloud), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_isolated_cutout.im".format(near_map),
		overwrite = True
			)    

	#imsubimage FAR map
	imsubimage(imagename="./{}_regrid.im".format(far_map), 	  
		mask="{}_cutout_regrid.im".format(cloud), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_isolated_cutout.im".format(far_map),
		overwrite = True
			)   
	"""





	"""

	###imsubimage FULL 7 to 9kpc of model to sailfish cutout and region###
	importfits(fitsimage="{}.fits".format(model_FULL),
			imagename="{}.im".format(model_FULL), 
			overwrite = True
				)  

	imsubimage(imagename="{}.im".format(model_FULL), 	  
		mask="{}_cutout_regrid.im".format(mask_name), 
		region = "{}_cutout_region.crtf".format(cloud),
		outfile="{}_model_7to9kpc_isolated_cutout.im".format(cloud),
		overwrite = True
			) 

	#export near/far masked clouds and full model of masked region

	"""






	"""
	exportfits(imagename="{}_isolated_cutout.im".format(near_map),
			fitsimage="{}_isolated_cutout.fits".format(near_map),
			overwrite=True
				)


	exportfits(imagename="{}_isolated_cutout.im".format(far_map),
			fitsimage="{}_isolated_cutout.fits".format(far_map),
			overwrite=True
				)




	exportfits(imagename="{}_model_7to9kpc_isolated_cutout.im".format(cloud),
			fitsimage="{}_model_7to9kpc_isolated_cutout.fits".format(cloud),
			overwrite=True
				)



	"""










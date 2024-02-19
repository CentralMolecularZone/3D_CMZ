"""
Uses CASA's imregrid task to regrid dendrogram mask to MALT09 and APEX CMZ
survey maps so that we can apply the masks later to analyse the spectral line
data within the dendrogram leaves.
"""
from astropy.io import fits

imregrid(imagename='./Masks/Final_continuum_clouds_mask.fits',
         template='./APEX_data/APEX_13CO_2014_merge.im',
         output='./Masks/Final_continuum_clouds_mask_regrid_to_APEX.im')

imregrid(imagename='./Masks/Final_continuum_clouds_mask.fits',
         template='./MALT90_data/CMZ_3mm_HNCO.im',
         output='./Masks/Final_continuum_clouds_mask_regrid_to_MALT90.im')

exportfits(imagename='./Masks/Final_continuum_clouds_mask_regrid_to_APEX.im',
           fitsimage='./Masks/Final_continuum_clouds_mask_regrid_to_APEX.fits')

exportfits(imagename='./Masks/Final_continuum_clouds_mask_regrid_to_MALT90.im',
           fitsimage='./Masks/Final_continuum_clouds_mask_regrid_to_MALT90.fits')


data = fits.getdata('./Masks/Final_continuum_clouds_mask_regrid_to_APEX.fits')
header = fits.getheader('./Masks/Final_continuum_clouds_mask_regrid_to_APEX.fits')
data[data>0] = 1
fits.writeto('./Masks/Final_continuum_clouds_mask_regrid_to_APEX_header_updated.fits',
            data, header)

data = fits.getdata('./Masks/Final_continuum_clouds_mask_regrid_to_MALT90.fits')
header = fits.getheader('./Masks/Final_continuum_clouds_mask_regrid_to_MALT90.fits')
data[data>0] = 1
fits.writeto('./Masks/Final_continuum_clouds_mask_regrid_to_MALT90_header_updated.fits',
            data, header)

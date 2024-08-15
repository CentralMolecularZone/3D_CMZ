from astropy.io import fits
import numpy as np
import glob
import os

file_list = [x for x in glob.glob('*_HNCO_cube_max.fits')]

for file in file_list:
    unique_ids = list(set([int(file.split('_')[0]) for file in file_list]))
    for id in unique_ids:
        if not os.path.exists(os.path.join('./weighted_masks/',os.path.splitext(id_files[i])[0] + '.mask.fits')):
            id_files = [file for file in file_list if int(file.split('_')[0]) == id]
            num_files = len(id_files)

            file_data = [fits.open(file)[0].data for file in id_files]
            headers = [fits.open(file)[0].header for file in id_files]

            sum_data = np.sum(file_data, axis=0)

            masks = []
            for i in range(num_files):
                weighted_mask = file_data[i] / sum_data
                masks.append(weighted_mask)

            for i, mask in enumerate(masks):
                mask_hdu = fits.PrimaryHDU(mask, headers[i])
                mask_hdulist = fits.HDUList([mask_hdu])
                mask_filename = os.path.splitext(id_files[i])[0] + '.mask.fits'
                mask_hdulist.writeto(os.path.join('./weighted_masks/', mask_filename))

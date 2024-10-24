import os
import numpy as np
import astrodendro
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
from astropy import wcs
import reproject
from paths import rootpath


lum_no70_fn = 'luminosity_gc_no70_conv36_with4pi-stupidunitsinheader.fits'
lum_fn = 'luminosity_warm_conv36-v8-nonans_beam_stupidunitsinheader.fits'
colfile = rootpath+'column_properunits_conv36_source_only.fits'

lum_fh = fits.open(lum_fn)
lum_no70_fh = fits.open(lum_no70_fn)
col_hdr = fits.getheader(colfile)

lum_wcs = wcs.WCS(lum_fh[0].header)
lum_no70_wcs = wcs.WCS(lum_no70_fh[0].header)
col_wcs = wcs.WCS(col_hdr)

lum_pix_area = wcs.utils.proj_plane_pixel_area(lum_wcs)
lum_no70_pix_area = wcs.utils.proj_plane_pixel_area(lum_no70_wcs)
col_pix_area = wcs.utils.proj_plane_pixel_area(col_wcs)

# convert units from Lsun/pix to Lsun/deg^2
lum_fh[0].data /= lum_pix_area
lum_no70_fh[0].data /= lum_no70_pix_area

# reproject both luminosity images to the column space
lum_data = reproject.reproject_interp(lum_fh, col_hdr)[0]
lum_no70_data = reproject.reproject_interp(lum_no70_fh, col_hdr)[0]

# convert back to Lsun units (because they're easier to work with; this could
# be done at any point)
lum_data *= col_pix_area
lum_no70_data *= col_pix_area
lum_data = u.Quantity(lum_data, u.L_sun)
lum_no70_data = u.Quantity(lum_no70_data, u.L_sun)

dendro_filename = 'cmz_column_dendrogram.hdf5'
dendro_filename = rootpath+"cmz_dendrogram_large.fits"
if os.path.exists(dendro_filename):
    dend = astrodendro.Dendrogram.load_from(dendro_filename)
else:
    raise ValueError("Run dendrogram_catalog first")

new_columns = {'peak_lum': [], 'mean_lum': [], 'total_lum': [], 'total_lum_SFR': [],
               'peak_lum_no70': [], 'mean_lum_no70': [], 'total_lum_no70': [], 'total_lum_no70_SFR': [],
              }

# Kennicutt relation from Ash's paper
lum_to_sfr_factor = 4.5e-44/(u.erg/u.s/u.M_sun*u.yr)

for struct_id in range(len(dend)):
    struct = dend[struct_id]
    mask = struct.get_mask()

    new_columns['peak_lum'].append(np.nanmax(lum_data[mask]))
    new_columns['mean_lum'].append(np.nanmean(lum_data[mask]))
    new_columns['total_lum'].append(np.nansum(lum_data[mask]))
    new_columns['total_lum_SFR'].append((np.nansum(lum_data[mask])*lum_to_sfr_factor).to(u.Msun/u.yr))
    new_columns['peak_lum_no70'].append(np.nanmax(lum_no70_data[mask]))
    new_columns['mean_lum_no70'].append(np.nanmean(lum_no70_data[mask]))
    new_columns['total_lum_no70'].append(np.nansum(lum_no70_data[mask]))
    new_columns['total_lum_no70_SFR'].append((np.nansum(lum_no70_data[mask])*lum_to_sfr_factor).to(u.Msun/u.yr))


cat = Table.read('column_and_temperature_withHCN_catalog.ipac',
                 format='ascii.ipac')

for key in new_columns:
    cat.add_column(Column(data=u.Quantity(new_columns[key]), name=key))

cat.write('column_and_temperature_and_luminosity_and_HCN_catalog.ipac',
          format='ascii.ipac', overwrite=True)

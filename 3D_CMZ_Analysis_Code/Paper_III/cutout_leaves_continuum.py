"""
Takes the dendrogram FITS file and cuts out sub-fields of the HiGAL continuum
data based on the leaf sizes.
"""

from __future__ import division
import numpy as np
import astrodendro
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from astropy.stats import mad_std
from astropy.table import Table, Column
from astropy.nddata import Cutout2D
plt.style.use('classic')
# _____________________________________

# Continuum file
file = "./../Data/higal_data/column_properunits_conv36_source_only.fits"
hdu = fits.open(file)[0]

# Some constants and stuff ...
pc          = const.pc
distance    = 8100*u.pc

# Some pixel area stuff ...
pix_width       = head['cdelt2']*u.deg
pix_width_pc    = (pix_width).to(u.rad).value * distance
pixels_1pc      = (1/pix_width_pc).value
pc_sc           = (pixels_1pc*pix_width).value

# Load dendrogram
dend = astrodendro.Dendrogram.load_from("./../Dendrogram_files/clouds_only_dendrogram.fits")
leaves = dend.leaves[9:(len(dend.leaves)-3)]

metadata = {}
metadata['data_unit'] = u.MJy / u.sr
metadata['spatial_scale'] = (np.abs(head['CDELT2'])) * u.deg
metadata['WCS'] = wcs

cat = pp_catalog(leaves, metadata)

# Cut out continuum based on radius from leaf area.
i=1
for leaf in leaves:
    hdu = fits.open(file)[0]
    wcs = WCS(hdu.header)
    position = (cat[cat['_idx']==leaf.idx]['x_cen'], cat[cat['_idx']==leaf.idx]['y_cen'])
    size = (4*int(np.around(np.sqrt(cat[cat['_idx']==leaf.idx]['area_exact']/np.pi)/pix_width.value)),
    4*int(np.around(np.sqrt(cat[cat['_idx']==leaf.idx]['area_exact']/np.pi)/pix_width.value)))
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)
    hdu.data = cutout.data
    hdu.header.update(cutout.wcs.to_header())
    hdu.writeto('./../Continuum_cutouts/'+str(i)+"_cutout.fits", overwrite=True)
    i=i+1

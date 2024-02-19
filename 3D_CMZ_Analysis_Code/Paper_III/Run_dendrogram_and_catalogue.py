"""
- Computes dendrogram fomr Herschel HiGAL dust column density map.
- Plots map with dendrogram leaves (i.e. clouds in this case) overlaid, with their
respective indicies.
- Outputs cloud catalogue + some quantities (e.g. dust temp, vlsr, delta_v) in
ascii & latex format.
"""

import os
import aplpy
import numpy as np
import astrodendro
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from astropy import constants as const
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic
from astropy import wcs
from astropy.wcs import WCS
from astropy.stats import mad_std
from astropy.table import Table, Column
import matplotlib.pyplot as plt
from astrodendro.analysis import PPStatistic, MetadataQuantity
plt.style.use('classic')

class MyPPStatistic(PPStatistic):
    """
    It turns out the existing PPStatistic class doesn't work very well with
    cm^-2 units, so I (credit Adam Ginsburg) made this new class to get the
    quantities we want.
    """

    distance = MetadataQuantity('distance', 'Distance to the target', strict=True)
    particle_mass = MetadataQuantity('particle_mass', 'Mass of each particle', strict=True)

    @property
    def average_column(self):

        average_col = self.stat.mom0() * self.data_unit / self.stat.count()

        return average_col.to(u.cm**-2)

    @property
    def mass(self):

        pixel_area = ((self.spatial_scale * self.distance)**2).to(u.cm**2,
        u.dimensionless_angles())
        mass = self.stat.mom0() * self.data_unit * pixel_area * self.particle_mass

        return mass.to(u.M_sun)

    @property
    def median_column(self):
        return np.nanmedian(self.stat.values) * self.data_unit

    @property
    def peak_column(self):
        return np.nanmax(self.stat.values) * self.data_unit
###

file = "./../Data/higal_data/column_properunits_conv36_source_only.fits"
head = fits.getheader(file)

# Some constants and stuff ...
pc          = const.pc
distance    = 8100*u.pc

# Some pixel area stuff ...
pix_width       = head['cdelt2']*u.deg
pix_width_pc    = (pix_width).to(u.rad).value * distance
pixels_1pc      = (1/pix_width_pc).value
pc_sc           = (pixels_1pc*pix_width).value

# Some file stuff ...
colfile = './../Data/higal_data/column_properunits_conv36_source_only.fits'
data = fits.getdata(colfile)
head = fits.getheader(colfile)
mywcs = wcs.WCS(fits.getheader(colfile))

# Load Brick dendrogram (had to be done separately to isolate the Brick)
dend_brick = astrodendro.Dendrogram.load_from("./../Dendrogram_files/separate_brick16_dendrogram.fits")
leaves_brick = dend_brick.leaves[-5]

# Compute dendrogram for HiGAL column density map
dend = astrodendro.Dendrogram.compute(data, wcs=mywcs, min_value=2e22, min_delta=5e22, min_npix=100)
dend.save_to('./../Dendrogram_files/clouds_only_dendrogram.hdf5')
dend.save_to('./../Dendrogram_files/clouds_only_dendrogram.fits')
leaves = dend.leaves[9:(len(dend.leaves)-3)]

# Get mask for dendrogram leaves ('clouds')
hdu = fits.open(file)[0]
mask = np.zeros(hdu.data.shape, dtype=bool)
for leaf in leaves:
    if leaf.idx == 45:
        mask = mask | leaves_brick.get_mask()
    else:
        mask = mask | leaf.get_mask()
mask_hdu = fits.PrimaryHDU(mask.astype('short'), hdu.header)

metadata = {}
metadata['data_unit'] = u.MJy / u.sr
metadata['spatial_scale'] = (np.abs(head['CDELT2'])) * u.deg
metadata['WCS'] = mywcs

# Force Brick catalogue data into full catalogue
cat = pp_catalog(leaves, metadata)
cat_brick = pp_catalog(dend_brick.leaves, metadata)
cat[cat['_idx']==45] = cat_brick[-5]
cat[cat['x_cen']==1032.6463805131666][0][0] = 15

# Plot HiGAL column density map with leave contours & labels
l = np.zeros(len(cat))
b = np.zeros(len(cat))
for i in range(0,len(cat),1):
    l[i] = cat[i][8]
    b[i] = cat[i][9]

fig=aplpy.FITSFigure(file)
lab = fig.pixel2world(l,b)
i = 0
for structure in leaves:
    fig.add_label(lab[0][i], lab[1][i], str(i+1), color='yellow',size=5)
    i = i+1
fig.recenter(0.45, -0.07, width=400*pc_sc, height=100*pc_sc)
fig.show_colorscale(cmap='Blues',pmin=0,pmax=99.9)
fig.add_scalebar(50*pc_sc)
fig.scalebar.show(50*pc_sc)
fig.scalebar.set_corner('top right')
fig.scalebar.set_color('black')
fig.scalebar.set_label('50 pc')
fig.scalebar.set_font(size=12)
fig.scalebar.set_linewidth(1)
fig.tick_labels.set_xformat('dd.dd')
fig.tick_labels.set_yformat('dd.dd')
fig.ticks.set_color('black')
fig.axis_labels.set_font(size=16)
fig.axis_labels.set_xpad(8)
fig.axis_labels.set_ypad(8)
fig.show_contour(mask_hdu, colors='red', linewidths=0.1)
fig.save('./../Figs/HiGAL_column_map_with_leaf_contours.eps')
plt.close()

# Build final catalogue
metadata = {}
metadata['data_unit'] = u.cm**-2
metadata['spatial_scale'] = mywcs.wcs.cdelt[1] * u.deg
metadata['beam_major'] = 36 * u.arcsec
metadata['beam_minor'] = 36 * u.arcsec
metadata['wcs'] = mywcs
metadata['distance'] = 8.1*u.kpc
metadata['particle_mass'] = 2.8*u.Da

fields = ['major_sigma', 'minor_sigma', 'radius', 'area_ellipse', 'area_exact',
          'position_angle', 'x_cen', 'y_cen', 'average_column',
          'median_column', 'peak_column', 'mass', ]


cat = astrodendro.analysis._make_catalog(structures=dend, fields=fields,
                                         metadata=metadata,
                                         statistic=MyPPStatistic,
                                         verbose=False)

cat_brick = astrodendro.analysis._make_catalog(structures=dend_brick.leaves,
                                               fields=fields,
                                               metadata=metadata,
                                               statistic=MyPPStatistic,
                                               verbose=False)

cat.rename_column('x_cen','l_cen')
cat.rename_column('y_cen','b_cen')

cat_leaves = cat[[structure.idx for structure in leaves]]
cat_leaves[cat_leaves['_idx']==45] = cat_brick[-5]
cat_leaves[cat_leaves['l_cen']==0.2546519427632299][0][0] = 15

temperature_filename = './../Data/higal_data/temp_conv36_source_only.fits'
temdata = fits.getdata(temperature_filename)

dendro_filename = './../Dendrogram_files/clouds_only_dendrogram.fits'
dend = astrodendro.Dendrogram.load_from(dendro_filename)

new_columns = {'peak_tem': [], 'mean_tem': [], 'median_tem': []}

# Append dust temperature data
for struct_id in range(len(leaves)):
    struct = leaves[struct_id]
    mask = struct.get_mask()

    new_columns['peak_tem'].append(np.nanmax(temdata[mask]))
    new_columns['mean_tem'].append(np.nanmean(temdata[mask]))
    new_columns['median_tem'].append(np.nanmedian(temdata[mask]))

for key in new_columns:
    cat_leaves.add_column(Column(data=new_columns[key]*u.K, name=key))

for i in range(len(cat_leaves)):
    cat_leaves['_idx'][i] = i+1

cat_leaves['radius'] = (cat_leaves['radius'].to(u.rad).value)*8100
cat_leaves['radius'].unit='pc'

cat_leaves['area_exact'] = (cat_leaves['area_exact'].to(u.rad*u.rad).value)*(8100**2)
cat_leaves['area_exact'].unit='pc$^{2}$'

# Compute equivalent radii
cat_leaves.add_column(Column(data=np.around(np.sqrt((cat_leaves['area_exact'] / np.pi) ),
                      decimals=1), name="Rad"))

general_props = ['_idx',
  'area_exact',
  'l_cen',
  'b_cen',
  'median_column',
  'peak_column',
  'mass',
  'Rad',
  'median_tem',
  'peak_tem']

cat_leaves_general = cat_leaves[general_props]

# Add vlsr, 2nd moment/linewidths, and colloquial cloud names.
v_map = {1: '19', 2: '39', 3: '16', 4: '-56', 5: '-29, -21', 6: '28, 58',
        7: '85', 8: '15', 9: '54', 10: '83', 11: '50', 12: '48', 13: '50',
        14: '62', 15: '35', 16: '51', 17: '49', 18: '29', 19: '53', 20: '21',
        21: '8, 39', 22: '-2'}

mom2_map = {1: '3', 2: '15', 4: '5', 5: '3', 6: '17', 7: '9', 8: '12', 9: '15',
            10: '20', 11: '8', 12: '10', 13: '5', 14: '12', 15: '11', 16: '6',
            17: '11', 18: '10', 19: '11', 20: '9', 21: '7', 22: '3'}

lw_map = {}

name_map = {14: 'Sagittarius B2', 13: 'G1.651-0.050', 17: 'G1.602+0.018',
            18: 'Dust ridge clouds E and F', 20: 'Dust ridge cloud D',
            21: 'Dust ridge cloud C', 22: 'Dust ridge cloud B', 15: 'The Brick',
            9: 'Straw and Sticks clouds', 11: 'Stone cloud',
            35: 'Three Little Pigs',12: '50 km/s cloud', 8: '20 km/s cloud',
            4: 'Sagittarius C'}

v_col = []
for row in cat_leaves_general:
    if row['_idx'] in v_map:
        v_col.append(v_map[row['_idx']])
    else:
        v_col.append('-')

mom2_col = []
for row in cat_leaves_general:
    if row['_idx'] in mom2_map:
        mom2_col.append(mom2_map[row['_idx']])
    else:
        mom2_col.append('-')

lw_col = []
for row in cat_leaves_general:
    if row['_idx'] in lw_map:
        lw_col.append(lw_map[row['_idx']])
    else:
        lw_col.append('-')

name_col = []
for row in cat_leaves_general:
    if row['_idx'] in name_map:
        name_col.append(name_map[row['_idx']])
    else:
        name_col.append('-')

cat_leaves_general.add_column(Table.Column(data=v_col, name='v_cen'))
cat_leaves_general.add_column(Table.Column(data=mom2_col, name='mean_mom2'))
cat_leaves_general.add_column(Table.Column(data=lw_col, name='fitted_lw'))
cat_leaves_general.add_column(Table.Column(data=name_col, name='Common_Name'))

# Format columns as desired.
cat_leaves_general['area_exact'].format = '6.0f'
cat_leaves_general['l_cen'].format = '6.3f'
cat_leaves_general['b_cen'].format = '6.3f'
cat_leaves_general['mass'].format = '%.1E'
cat_leaves_general['median_column'].format = '%.1E'
cat_leaves_general['peak_column'].format = '%.1E'
cat_leaves_general['median_tem'].format = '6.0f'
cat_leaves_general['peak_tem'].format = '6.0f'

# Output catalogue
cat_leaves_general.write('cloud_only_catalog_with_temp.ipac',
                         format='ascii.ipac', overwrite=True)
cat_leaves_general.write('cloud_only_catalog_with_temp.tex',
                         format='ascii.latex', overwrite=True)

import os
import numpy as np
import astrodendro
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
from paths import rootpath

temperature_filename = rootpath+'temp_conv36_source_only.fits'

temdata = fits.getdata(temperature_filename)

dendro_filename = 'cmz_column_dendrogram.hdf5'
dendro_filename = rootpath+"cmz_dendrogram_large.fits"
if os.path.exists(dendro_filename):
    dend = astrodendro.Dendrogram.load_from(dendro_filename)
else:
    raise ValueError("Run dendrogram_catalog first")

assert dend.data.shape == temdata.shape,"Temperature data should be same shape as dendrogram data"

new_columns = {'peak_tem': [], 'mean_tem': [], 'median_tem': []}

for struct_id in range(len(dend)):
    struct = dend[struct_id]
    mask = struct.get_mask()

    new_columns['peak_tem'].append(np.nanmax(temdata[mask]))
    new_columns['mean_tem'].append(np.nanmean(temdata[mask]))
    new_columns['median_tem'].append(np.nanmedian(temdata[mask]))


cat = Table.read('column_catalog.ipac', format='ascii.ipac')

for key in new_columns:
    cat.add_column(Column(data=new_columns[key]*u.K, name=key))

cat.write('column_and_temperature_catalog.ipac', format='ascii.ipac', overwrite=True)

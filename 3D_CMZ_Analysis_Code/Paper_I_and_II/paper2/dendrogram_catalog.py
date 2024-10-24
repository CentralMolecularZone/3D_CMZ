import os
import numpy as np
import astrodendro
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astrodendro.analysis import PPStatistic, MetadataQuantity
#from paths import rootpath

rootpath='/Users/battersby/Dropbox/higal_cmz_dendro/AA_re_extraction_test_July2022/'

colfile = rootpath+'column_properunits_conv36_source_only.fits'
data = fits.getdata(colfile)
mywcs = wcs.WCS(fits.getheader(colfile))

dendro_filename = 'cmz_column_dendrogram.hdf5'
dendro_filename = rootpath+"cmz_dendrogram_large.fits"
if os.path.exists(dendro_filename):
    dend = astrodendro.Dendrogram.load_from(dendro_filename)
else:
    dend = astrodendro.Dendrogram.compute(data, wcs=mywcs, min_value=2e22,
                                          min_delta=5e22, min_npix=10)

    dend.save_to(dendro_filename)

    maskdata = np.zeros_like(dend.data)
    for structure in dend.all_structures:
        maskdata[structure.get_mask()] = structure.idx
    fh = fits.open(colfile)
    fh[0].data = maskdata.astype('int')
    fh.writeto('cmz_column_dendrogram_mask.fits', overwrite=True)

metadata = {}
metadata['data_unit'] = u.cm**-2
metadata['spatial_scale'] = mywcs.wcs.cdelt[1] * u.deg
metadata['beam_major'] = 36 * u.arcsec
metadata['beam_minor'] = 36 * u.arcsec
metadata['wcs'] = mywcs
metadata['distance'] = 8.1*u.kpc #updated July 2022
metadata['particle_mass'] = 2.8*u.Da # 2.8 AMU (Dalton=AMU), per Kauffmann 2008 appendix

fields = ['major_sigma', 'minor_sigma', 'radius', 'area_ellipse', 'area_exact',
          'position_angle', 'x_cen', 'y_cen', 'average_column',
          'median_column', 'peak_column', 'mass', ]


class MyPPStatistic(PPStatistic):
    """
    It turns out the existing PPStatistic class doesn't work very well with
    cm^-2 units, so I made this new class to get the quantities we want.
    """

    distance = MetadataQuantity('distance', 'Distance to the target', strict=True)
    particle_mass = MetadataQuantity('particle_mass', 'Mass of each particle', strict=True)

    @property
    def average_column(self):

        # mom0 sums the column pixels, data unit should be surface density,
        # stat.count is # of pixels
        average_col = self.stat.mom0() * self.data_unit / self.stat.count()

        return average_col.to(u.cm**-2)

    @property
    def mass(self):

        pixel_area = ((self.spatial_scale * self.distance)**2).to(u.cm**2, u.dimensionless_angles())
        mass = self.stat.mom0() * self.data_unit * pixel_area * self.particle_mass

        return mass.to(u.M_sun)

    @property
    def median_column(self):
        return np.nanmedian(self.stat.values) * self.data_unit

    @property
    def peak_column(self):
        return np.nanmax(self.stat.values) * self.data_unit


#cat = astrodendro.pp_catalog(dend, metadata=metadata, fields=fields,
#                             statistic=MyPPStatistic)
cat = astrodendro.analysis._make_catalog(structures=dend, fields=fields,
                                         metadata=metadata,
                                         statistic=MyPPStatistic,
                                         verbose=False)
cat.rename_column('x_cen','l_cen')
cat.rename_column('y_cen','b_cen')

cat.write('column_catalog.ipac', format='ascii.ipac', overwrite=True)

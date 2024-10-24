import os
import pyspeckit
import numpy as np
import astrodendro
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy.table import Table, Column
from astropy.stats import mad_std
from spectral_cube import SpectralCube
import spectral_cube
from paths import rootpath
from astropy import log
#import paths

SIGMA2FWHM = 2. * np.sqrt(2. * np.log(2.))

dendro_filename = 'cmz_column_dendrogram.hdf5'
dendro_filename = rootpath+"cmz_dendrogram_large.fits"
if os.path.exists(dendro_filename):
    dend = astrodendro.Dendrogram.load_from(dendro_filename)
else:
    raise ValueError("Run dendrogram_catalog first")

extract_spectra = False
# check whether spectral extraction is all done
for structure in dend.all_structures:
    hnco_specname = 'spectra/HNCO_spectrum_idx{0}.fits'.format(structure.idx)
    hc3n_specname = 'spectra/HC3N_spectrum_idx{0}.fits'.format(structure.idx)
    hcn_specname = 'spectra/HCN_spectrum_idx{0}.fits'.format(structure.idx)
    if not os.path.exists(hcn_specname):
        print("No archived {0}".format(hcn_specname))
        extract_spectra = True
    if not os.path.exists(hc3n_specname):
        print("No archived {0}".format(hc3n_specname))
        extract_spectra = True
    if not os.path.exists(hnco_specname):
        print("No archived {0}".format(hnco_specname))
        extract_spectra = True

extract_spectra = True
#extract_spectra = False # set manually b/c searching for files isn't good enough
# when NaN regions are skipped
if extract_spectra and not('cube' in locals() and 'cube_hc3n' in locals() and 'cube_hnco' in locals()):
    cube = SpectralCube.read('CMZ_3mm_HCN.fits')
    med_hcn = cube.with_mask(cube < 0.1*cube.unit).median(axis=0)
    cube = cube - med_hcn

    cube_hc3n = SpectralCube.read('CMZ_3mm_HC3N.fits')
    med_hc3n = cube_hc3n.with_mask(cube_hc3n < 0.1*cube_hc3n.unit).median(axis=0)
    cube_hc3n = cube_hc3n - med_hc3n

    cube_hnco = SpectralCube.read('CMZ_3mm_HNCO.fits')
    med_hnco = cube_hnco.with_mask(cube_hnco < 0.1*cube_hnco.unit).median(axis=0)
    cube_hnco = cube_hnco - med_hnco

    cube_twelveco = cube_12co = SpectralCube.read('gc_12co.fits')
    cube_thirteenco = cube_13co = SpectralCube.read('gc_13co.fits')

    colfile = rootpath+'column_properunits_conv36_source_only.fits'
    header = fits.getheader(colfile)
    mywcs = wcs.WCS(header)


    cube_header = cube.header.copy()
    # in theory, this should update only the spatial component of the header
    cube_header.update(mywcs.to_header())
    cube_header['NAXIS1'] = header['NAXIS1']
    cube_header['NAXIS2'] = header['NAXIS2']

    # This is computationally VERY expensive...
    # have to do some serious hacks to get it down to a feasible size...
    #reproj_cube = cube.reproject(cube_header)

    # these are the outer extents that overlap with the HCN image
    lx,ly,ux,uy = 805,466,1978,720
    lx,ly,ux,uy = 544,320,1342,489
    naxis1 = ux-lx
    naxis2 = uy-ly
    crpix1 = header['CRPIX1'] - lx
    crpix2 = header['CRPIX2'] - ly
    cube_header['NAXIS1'] = naxis1
    cube_header['NAXIS2'] = naxis2
    cube_header['CRPIX1'] = crpix1
    cube_header['CRPIX2'] = crpix2

    reproj_cube = cube.reproject(cube_header)

    cube_header_hc3n = cube_hc3n.header.copy()
    # in theory, this should update only the spatial component of the header
    cube_header_hc3n.update(mywcs.to_header())
    cube_header_hc3n['NAXIS1'] = header['NAXIS1']
    cube_header_hc3n['NAXIS2'] = header['NAXIS2']

    cube_header_hc3n['NAXIS1'] = naxis1
    cube_header_hc3n['NAXIS2'] = naxis2
    cube_header_hc3n['CRPIX1'] = crpix1
    cube_header_hc3n['CRPIX2'] = crpix2

    reproj_cube_hc3n = cube_hc3n.reproject(cube_header_hc3n)

    cube_header_hnco = cube_hnco.header.copy()
    # in theory, this should update only the spatial component of the header
    cube_header_hnco.update(mywcs.to_header())
    cube_header_hnco['NAXIS1'] = header['NAXIS1']
    cube_header_hnco['NAXIS2'] = header['NAXIS2']

    cube_header_hnco['NAXIS1'] = naxis1
    cube_header_hnco['NAXIS2'] = naxis2
    cube_header_hnco['CRPIX1'] = crpix1
    cube_header_hnco['CRPIX2'] = crpix2

    reproj_cube_hnco = cube_hnco.reproject(cube_header_hnco)

    cube_header_co = cube_header_hnco.copy()
    cube_header_co['CTYPE3'] = cube_12co.header['CTYPE3']
    cube_header_co['NAXIS3'] = cube_12co.header['NAXIS3']
    reproj_cube_twelveco = cube_12co.reproject(cube_header_co)
    cube_header_co['CTYPE1'] = cube_13co.header['CTYPE1']
    cube_header_co['CTYPE2'] = cube_13co.header['CTYPE2']
    cube_header_co['CTYPE3'] = cube_13co.header['CTYPE3']
    cube_header_co['NAXIS3'] = cube_13co.header['NAXIS3']
    reproj_cube_thirteenco = cube_13co.reproject(cube_header_co)


new_columns = {'HCN_mom0': [], 'HCN_mom1': [], 'HCN_fwhm':[], 'HCN_spectrum': [],
               'HC3N_mom0': [], 'HC3N_mom1': [],'HC3N_fwhm':[],  'HC3N_spectrum': [],
               'HNCO_mom0': [], 'HNCO_mom1': [],'HNCO_fwhm':[],  'HNCO_spectrum': [],
               '12CO_mom0': [], '12CO_mom1': [],'12CO_fwhm':[],  '12CO_spectrum': [],
               '13CO_mom0': [], '13CO_mom1': [],'13CO_fwhm':[],  '13CO_spectrum': [],
               'has_zone_outside_data': [],
               'npix_total': [],
               'npix_cubes': [],
              }
units = {'HCN_mom0': u.K*u.km/u.s, 'HCN_mom1': u.km/u.s, 'HCN_fwhm': u.km/u.s,
         'HC3N_mom0': u.K*u.km/u.s, 'HC3N_mom1': u.km/u.s, 'HC3N_fwhm': u.km/u.s,
         'HNCO_mom0': u.K*u.km/u.s, 'HNCO_mom1': u.km/u.s, 'HNCO_fwhm': u.km/u.s,
         '12CO_mom0': u.K*u.km/u.s, '12CO_mom1': u.km/u.s, '12CO_fwhm': u.km/u.s,
         '13CO_mom0': u.K*u.km/u.s, '13CO_mom1': u.km/u.s, '13CO_fwhm': u.km/u.s,
         'has_zone_outside_data': None,
         'npix_total': None,
         'npix_cubes': None,
        }

for structure_indx in range(len(dend)):
    structure = dend[structure_indx]
    dend_obj_mask = structure.get_mask()
    npix = dend_obj_mask.sum()

    # hack to speed things up a LOT: crop more.
    dend_inds = structure.indices()

    view = [slice(dend_inds[0].min(), dend_inds[0].max()+1),
            slice(dend_inds[1].min(), dend_inds[1].max()+1)]
    submask = dend_obj_mask[view]

    has_zone_outside_data = False

    hcn_specname = 'spectra/HCN_spectrum_idx{0}.fits'.format(structure.idx)
    hc3n_specname = 'spectra/HC3N_spectrum_idx{0}.fits'.format(structure.idx)
    hnco_specname = 'spectra/HNCO_spectrum_idx{0}.fits'.format(structure.idx)
    twelveco_specname = 'spectra/12CO_spectrum_idx{0}.fits'.format(structure.idx)
    thirteenco_specname = 'spectra/13CO_spectrum_idx{0}.fits'.format(structure.idx)
    if extract_spectra:
        #((not os.path.exists(hcn_specname) or not os.path.exists(hc3n_specname) or not os.path.exists(hnco_specname)) and extract_spectra):

        # crop the cube differently from the mask because the cube is pre-cropped
        cubeview = [slice(None),
                    slice(dend_inds[0].min()-ly, dend_inds[0].max()+1-ly),
                    slice(dend_inds[1].min()-lx, dend_inds[1].max()+1-lx)]
        if (((cubeview[1].start < 0 and cubeview[1].stop < 0) or (cubeview[2].start < 0 and cubeview[2].stop < 0) or
             (cubeview[1].start > reproj_cube.shape[1]) or (cubeview[2].start > reproj_cube.shape[2]))):
            print("Skipped extracting source {0} because it was out of range".format(structure.idx))
            has_zone_outside_data = True
            new_columns['HCN_mom0'].append(np.nan)
            new_columns['HCN_mom1'].append(np.nan)
            new_columns['HCN_fwhm'].append(np.nan)
            new_columns['HCN_spectrum'].append(None)
            new_columns['HC3N_mom0'].append(np.nan)
            new_columns['HC3N_mom1'].append(np.nan)
            new_columns['HC3N_fwhm'].append(np.nan)
            new_columns['HC3N_spectrum'].append(None)
            new_columns['HNCO_mom0'].append(np.nan)
            new_columns['HNCO_mom1'].append(np.nan)
            new_columns['HNCO_fwhm'].append(np.nan)
            new_columns['HNCO_spectrum'].append(None)
            new_columns['12CO_mom0'].append(np.nan)
            new_columns['12CO_mom1'].append(np.nan)
            new_columns['12CO_fwhm'].append(np.nan)
            new_columns['12CO_spectrum'].append(None)
            new_columns['13CO_mom0'].append(np.nan)
            new_columns['13CO_mom1'].append(np.nan)
            new_columns['13CO_fwhm'].append(np.nan)
            new_columns['13CO_spectrum'].append(None)
            new_columns['has_zone_outside_data'].append(has_zone_outside_data)
            new_columns['npix_total'].append(npix)
            new_columns['npix_cubes'].append(0)
            print(cubeview, reproj_cube.shape)
            if structure.idx == 9:
                raise ValueError
            print("Skipped extracting source {0} because it was out of range".format(structure.idx))
            continue

        print("Structure {0}".format(structure.idx))
        print("old view: {0}".format(view))
        print("old cubeview: {0}".format(cubeview))
        if cubeview[1].start < 0:
            view[0] = slice(view[0].start-cubeview[1].start, view[0].stop)
            cubeview[1] = slice(0, cubeview[1].stop)
            has_zone_outside_data = True
        if cubeview[2].start < 0:
            view[1] = slice(view[1].start-cubeview[2].start, view[1].stop)
            cubeview[2] = slice(0, cubeview[2].stop)
            has_zone_outside_data = True
        if cubeview[2].stop > reproj_cube.shape[2]:
            view[1] = slice(view[1].start, view[1].stop - (cubeview[2].stop - reproj_cube.shape[2]))
            cubeview[2] = slice(cubeview[2].start, None)
            has_zone_outside_data = True
        if cubeview[1].stop > reproj_cube.shape[1]:
            view[0] = slice(view[0].start, view[0].stop - (cubeview[1].stop - reproj_cube.shape[1]))
            cubeview[1] = slice(cubeview[1].start, None)
            has_zone_outside_data = True
        print("new view: {0}".format(view))
        print("new cubeview: {0}".format(cubeview))
        print('original submask shape: {0}'.format(submask.shape))
        submask = dend_obj_mask[view]
        print('new submask shape: {0}'.format(submask.shape))
        assert not any(x==0 for x in submask.shape)
        print('reproj_cube shape: {0}'.format(reproj_cube[cubeview].shape))

    elif (not os.path.exists(hcn_specname) or not os.path.exists(hc3n_specname) or not os.path.exists(hnco_specname)):
        new_columns['HCN_mom0'].append(np.nan)
        new_columns['HCN_mom1'].append(np.nan)
        new_columns['HCN_fwhm'].append(np.nan)
        new_columns['HCN_spectrum'].append(None)
        new_columns['HC3N_mom0'].append(np.nan)
        new_columns['HC3N_mom1'].append(np.nan)
        new_columns['HC3N_fwhm'].append(np.nan)
        new_columns['HC3N_spectrum'].append(None)
        new_columns['HNCO_mom0'].append(np.nan)
        new_columns['HNCO_mom1'].append(np.nan)
        new_columns['HNCO_fwhm'].append(np.nan)
        new_columns['HNCO_spectrum'].append(None)
        new_columns['12CO_mom0'].append(np.nan)
        new_columns['12CO_mom1'].append(np.nan)
        new_columns['12CO_fwhm'].append(np.nan)
        new_columns['12CO_spectrum'].append(None)
        new_columns['13CO_mom0'].append(np.nan)
        new_columns['13CO_mom1'].append(np.nan)
        new_columns['13CO_fwhm'].append(np.nan)
        new_columns['13CO_spectrum'].append(None)
        print("Skipped source {0} because it has no file and extract=False".format(structure.idx))
        continue

    if os.path.exists(hcn_specname):
        header = fits.getheader(hcn_specname)
        meanspec = spectral_cube.lower_dimensional_structures.OneDSpectrum(
            value=fits.getdata(hcn_specname), wcs=wcs.WCS(header),
            unit=header['BUNIT'])
        dx_hcn = np.abs(fits.getheader(hcn_specname)['CDELT1'])*u.m/u.s
    else:
        cropcube = reproj_cube[cubeview].with_mask(submask[None,:,:])
        meanspec = cropcube.mean(axis=(1,2))
        assert meanspec.size == cube.shape[0]
        meanspec.write(hcn_specname, overwrite=True)
        print("Wrote spectrum to {0} with max {1}".format(hcn_specname, np.nanmax(meanspec)))
        dx_hcn = cropcube._pix_size_slice(0) * cropcube.spectral_axis.unit

    if os.path.exists(hc3n_specname):
        header = fits.getheader(hc3n_specname)
        meanspec_hc3n = spectral_cube.lower_dimensional_structures.OneDSpectrum(
            value=fits.getdata(hc3n_specname), wcs=wcs.WCS(header),
            unit=header['BUNIT'])
        dx_hc3n = np.abs(fits.getheader(hc3n_specname)['CDELT1'])*u.m/u.s
    else:
        cropcube_hc3n = reproj_cube_hc3n[cubeview].with_mask(submask[None,:,:])
        meanspec_hc3n = cropcube_hc3n.mean(axis=(1,2))
        assert meanspec_hc3n.size == cube_hc3n.shape[0]
        meanspec_hc3n.write(hc3n_specname, overwrite=True)
        print("Wrote spectrum to {0} with max {1}".format(hc3n_specname, np.nanmax(meanspec_hc3n)))
        dx_hc3n = cropcube_hc3n._pix_size_slice(0) * cropcube_hc3n.spectral_axis.unit

    if os.path.exists(hnco_specname):
        header = fits.getheader(hnco_specname)
        meanspec_hnco = spectral_cube.lower_dimensional_structures.OneDSpectrum(
            value=fits.getdata(hnco_specname), wcs=wcs.WCS(header),
            unit=header['BUNIT'])
        dx_hnco = np.abs(fits.getheader(hnco_specname)['CDELT1'])*u.m/u.s
    else:
        cropcube_hnco = reproj_cube_hnco[cubeview].with_mask(submask[None,:,:])
        meanspec_hnco = cropcube_hnco.mean(axis=(1,2))
        assert meanspec_hnco.size == cube_hnco.shape[0]
        meanspec_hnco.write(hnco_specname, overwrite=True)
        print("Wrote spectrum to {0} with max {1}".format(hnco_specname, np.nanmax(meanspec_hnco)))
        dx_hnco = cropcube_hnco._pix_size_slice(0) * cropcube_hnco.spectral_axis.unit

    if os.path.exists(twelveco_specname):
        header = fits.getheader(twelveco_specname)
        meanspec_twelveco = spectral_cube.lower_dimensional_structures.OneDSpectrum(
            value=fits.getdata(twelveco_specname), wcs=wcs.WCS(header),
            unit=header['BUNIT'])
        dx_twelveco = np.abs(fits.getheader(twelveco_specname)['CDELT1'])*u.m/u.s
    else:
        cropcube_twelveco = reproj_cube_twelveco[cubeview].with_mask(submask[None,:,:])
        meanspec_twelveco = cropcube_twelveco.mean(axis=(1,2))
        assert meanspec_twelveco.size == cube_twelveco.shape[0]
        meanspec_twelveco.write(twelveco_specname, overwrite=True)
        print("Wrote spectrum to {0} with max {1}".format(twelveco_specname, np.nanmax(meanspec_twelveco)))
        dx_twelveco = cropcube_twelveco._pix_size_slice(0) * cropcube_twelveco.spectral_axis.unit


    if os.path.exists(thirteenco_specname):
        header = fits.getheader(thirteenco_specname)
        meanspec_thirteenco = spectral_cube.lower_dimensional_structures.OneDSpectrum(
            value=fits.getdata(thirteenco_specname), wcs=wcs.WCS(header),
            unit=header['BUNIT'])
        dx_thirteenco = np.abs(fits.getheader(thirteenco_specname)['CDELT1'])*u.m/u.s
    else:
        cropcube_thirteenco = reproj_cube_thirteenco[cubeview].with_mask(submask[None,:,:])
        meanspec_thirteenco = cropcube_thirteenco.mean(axis=(1,2))
        assert meanspec_thirteenco.size == cube_thirteenco.shape[0]
        meanspec_thirteenco.write(thirteenco_specname, overwrite=True)
        print("Wrote spectrum to {0} with max {1}".format(thirteenco_specname, np.nanmax(meanspec_thirteenco)))
        dx_thirteenco = cropcube_thirteenco._pix_size_slice(0) * cropcube_thirteenco.spectral_axis.unit

    # only look within +/- max_width km/s of peak
    max_width = 80*u.km/u.s


    argmax = meanspec.argmax()
    vmax = meanspec.spectral_axis[argmax]
    mask = (meanspec > mad_std(meanspec)) & (meanspec.spectral_axis-vmax>-max_width) & (meanspec.spectral_axis-vmax < max_width)
    med = np.nanmedian(meanspec[~mask])
    meanspec_bl = meanspec - med
    hcn_mom0 = u.Quantity(meanspec_bl[mask].sum() * dx_hcn)
    hcn_mom1 = u.Quantity((meanspec_bl[mask] * meanspec.spectral_axis[mask]).sum() / meanspec_bl[mask].sum())
    hcn_mom2 = u.Quantity((meanspec_bl[mask] * (meanspec.spectral_axis[mask]-hcn_mom1)**2).sum() / meanspec_bl[mask].sum())**0.5 * SIGMA2FWHM

    sp_hcn = pyspeckit.Spectrum.from_hdu(meanspec.hdu)
    sp_hcn.xarr.convert_to_unit(u.km/u.s)
    sp_hcn.specfit(guesses=[meanspec[argmax].value, vmax.to(u.km/u.s).value, hcn_mom2.to(u.km/u.s).value])

    new_columns['HCN_mom0'].append(hcn_mom0.to(u.K*u.km/u.s).value)
    #new_columns['HCN_mom1'].append(hcn_mom1.to(u.km/u.s).value)
    #new_columns['HCN_fwhm'].append(hcn_mom2.to(u.km/u.s).value)
    new_columns['HCN_mom1'].append(sp_hcn.specfit.parinfo.SHIFT0.value)
    new_columns['HCN_fwhm'].append(sp_hcn.specfit.parinfo.WIDTH0.value * SIGMA2FWHM)
    new_columns['HCN_spectrum'].append(meanspec)

    hc3nargmax = meanspec_hc3n.argmax()
    hc3nvmax = meanspec_hc3n.spectral_axis[hc3nargmax]
    hc3nmask = (meanspec_hc3n > mad_std(meanspec_hc3n)) & (meanspec_hc3n.spectral_axis-hc3nvmax>-max_width) & (meanspec_hc3n.spectral_axis-hc3nvmax < max_width)
    #hc3nmask = (meanspec_hc3n > mad_std(meanspec_hc3n)) & (meanspec_hc3n.spectral_axis>-200*u.km/u.s) & (meanspec_hc3n.spectral_axis < 200 *u.km/u.s)
    hc3nmed = np.nanmedian(meanspec_hc3n[~hc3nmask])
    meanspec_hc3n_bl = meanspec_hc3n - hc3nmed
    hc3n_mom0 = u.Quantity(meanspec_hc3n_bl[hc3nmask].sum() * dx_hc3n)
    hc3n_mom1 = u.Quantity((meanspec_hc3n_bl[hc3nmask] * meanspec_hc3n.spectral_axis[hc3nmask]).sum() / meanspec_hc3n_bl[hc3nmask].sum())
    hc3n_mom2 = u.Quantity((meanspec_hc3n_bl[hc3nmask] * (meanspec_hc3n.spectral_axis[hc3nmask]-hc3n_mom1)**2).sum() / meanspec_hc3n_bl[hc3nmask].sum())**0.5 * SIGMA2FWHM
    if np.abs(hc3n_mom1 - hc3nvmax) > 30 * u.km/u.s:
        log.critical("Major disagreement between peak and moment 1.  Recomputing with a much narrower range.")
        hc3nmask = (meanspec_hc3n > mad_std(meanspec_hc3n)) & (meanspec_hc3n.spectral_axis-hc3nvmax>-30*u.km/u.s) & (meanspec_hc3n.spectral_axis-hc3nvmax < 30*u.km/u.s)
        hc3nmed = np.nanmedian(meanspec_hc3n[~hc3nmask])
        meanspec_hc3n_bl = meanspec_hc3n - hc3nmed
        hc3n_mom0 = u.Quantity(meanspec_hc3n_bl[hc3nmask].sum() * dx_hc3n)
        hc3n_mom1 = u.Quantity((meanspec_hc3n_bl[hc3nmask] * meanspec_hc3n.spectral_axis[hc3nmask]).sum() / meanspec_hc3n_bl[hc3nmask].sum())
        hc3n_mom2 = u.Quantity((meanspec_hc3n_bl[hc3nmask] * (meanspec_hc3n.spectral_axis[hc3nmask]-hc3n_mom1)**2).sum() / meanspec_hc3n_bl[hc3nmask].sum())**0.5 * SIGMA2FWHM

    #If we want to fit the brightest line, we can use this
    sp_hc3n = pyspeckit.Spectrum.from_hdu(meanspec_hc3n.hdu)
    sp_hc3n.xarr.convert_to_unit(u.km/u.s)
    #sp.plotter(figure=pl.figure(1))
    sp_hc3n.specfit(guesses=[meanspec_hc3n[hc3nargmax].value, hc3nvmax.to(u.km/u.s).value, hc3n_mom2.to(u.km/u.s).value])

    new_columns['HC3N_mom0'].append(hc3n_mom0.to(u.K*u.km/u.s).value)
    #new_columns['HC3N_mom1'].append(hc3n_mom1.to(u.km/u.s).value)
    #new_columns['HC3N_fwhm'].append(hc3n_mom2.to(u.km/u.s).value)
    new_columns['HC3N_mom1'].append(sp_hc3n.specfit.parinfo.SHIFT0.value)
    new_columns['HC3N_fwhm'].append(sp_hc3n.specfit.parinfo.WIDTH0.value * SIGMA2FWHM)
    new_columns['HC3N_spectrum'].append(meanspec_hc3n)

    hncoargmax = meanspec_hnco.argmax()
    hncovmax = meanspec_hnco.spectral_axis[hncoargmax]
    hncomask = (meanspec_hnco > mad_std(meanspec_hnco)) & (meanspec_hnco.spectral_axis-hncovmax>-max_width) & (meanspec_hnco.spectral_axis-hncovmax < max_width)
    #hncomask = (meanspec_hnco > mad_std(meanspec_hnco)) & (meanspec_hnco.spectral_axis>-200*u.km/u.s) & (meanspec_hnco.spectral_axis < 200 *u.km/u.s)
    hncomed = np.nanmedian(meanspec_hnco[~hncomask])
    meanspec_hnco_bl = meanspec_hnco - hncomed
    hnco_mom0 = u.Quantity(meanspec_hnco_bl[hncomask].sum() * dx_hnco)
    hnco_mom1 = u.Quantity((meanspec_hnco_bl[hncomask] * meanspec_hnco.spectral_axis[hncomask]).sum() / meanspec_hnco_bl[hncomask].sum())
    hnco_mom2 = u.Quantity((meanspec_hnco_bl[hncomask] * (meanspec_hnco.spectral_axis[hncomask]-hnco_mom1)**2).sum() / meanspec_hnco_bl[hncomask].sum())**0.5 * SIGMA2FWHM

    #If we want to fit the brightest line, we can use this
    sp_hnco = pyspeckit.Spectrum.from_hdu(meanspec_hnco.hdu)
    sp_hnco.xarr.convert_to_unit(u.km/u.s)
    #sp.plotter(figure=pl.figure(1))
    sp_hnco.specfit(guesses=[meanspec_hnco[hncoargmax].value, vmax.to(u.km/u.s).value, hnco_mom2.to(u.km/u.s).value])

    new_columns['HNCO_mom0'].append(hnco_mom0.to(u.K*u.km/u.s).value)
    #new_columns['HNCO_mom1'].append(hnco_mom1.to(u.km/u.s).value)
    #new_columns['HNCO_fwhm'].append(hnco_mom2.to(u.km/u.s).value)
    new_columns['HNCO_mom1'].append(sp_hnco.specfit.parinfo.SHIFT0.value)
    new_columns['HNCO_fwhm'].append(sp_hnco.specfit.parinfo.WIDTH0.value * SIGMA2FWHM)
    new_columns['HNCO_spectrum'].append(meanspec_hnco)

    twelvecoargmax = meanspec_twelveco.argmax()
    if np.isfinite(meanspec_twelveco.max()):
        twelvecovmax = meanspec_twelveco.spectral_axis[twelvecoargmax]
        twelvecomask = (meanspec_twelveco > mad_std(meanspec_twelveco)) & (meanspec_twelveco.spectral_axis-twelvecovmax>-max_width) & (meanspec_twelveco.spectral_axis-twelvecovmax < max_width)
        #twelvecomask = (meanspec_twelveco > mad_std(meanspec_twelveco)) & (meanspec_twelveco.spectral_axis>-200*u.km/u.s) & (meanspec_twelveco.spectral_axis < 200 *u.km/u.s)
        twelvecomed = np.nanmedian(meanspec_twelveco[~twelvecomask])
        meanspec_twelveco_bl = meanspec_twelveco - twelvecomed
        twelveco_mom0 = u.Quantity(meanspec_twelveco_bl[twelvecomask].sum() * dx_twelveco)
        twelveco_mom1 = u.Quantity((meanspec_twelveco_bl[twelvecomask] * meanspec_twelveco.spectral_axis[twelvecomask]).sum() / meanspec_twelveco_bl[twelvecomask].sum())
        twelveco_mom2 = u.Quantity((meanspec_twelveco_bl[twelvecomask] * (meanspec_twelveco.spectral_axis[twelvecomask]-twelveco_mom1)**2).sum() / meanspec_twelveco_bl[twelvecomask].sum())**0.5 * SIGMA2FWHM

        #If we want to fit the brightest line, we can use this
        sp_twelveco = pyspeckit.Spectrum.from_hdu(meanspec_twelveco.hdu)
        sp_twelveco.xarr.convert_to_unit(u.km/u.s)
        #sp.plotter(figure=pl.figure(1))
        sp_twelveco.specfit(guesses=[meanspec_twelveco[twelvecoargmax].value, vmax.to(u.km/u.s).value, twelveco_mom2.to(u.km/u.s).value])

        new_columns['12CO_mom0'].append(twelveco_mom0.to(u.K*u.km/u.s).value)
        #new_columns['12CO_mom1'].append(twelveco_mom1.to(u.km/u.s).value)
        #new_columns['12CO_fwhm'].append(twelveco_mom2.to(u.km/u.s).value)
        new_columns['12CO_mom1'].append(sp_twelveco.specfit.parinfo.SHIFT0.value)
        new_columns['12CO_fwhm'].append(sp_twelveco.specfit.parinfo.WIDTH0.value * SIGMA2FWHM)
        new_columns['12CO_spectrum'].append(meanspec_twelveco)
    else:
        new_columns['12CO_mom0'].append(np.nan)
        new_columns['12CO_mom1'].append(np.nan)
        new_columns['12CO_fwhm'].append(np.nan)
        new_columns['12CO_spectrum'].append(None)

    thirteencoargmax = meanspec_thirteenco.argmax()
    if np.isfinite(meanspec_thirteenco.max()):
        thirteencovmax = meanspec_thirteenco.spectral_axis[thirteencoargmax]
        thirteencomask = (meanspec_thirteenco > mad_std(meanspec_thirteenco)) & (meanspec_thirteenco.spectral_axis-thirteencovmax>-max_width) & (meanspec_thirteenco.spectral_axis-thirteencovmax < max_width)
        #thirteencomask = (meanspec_thirteenco > mad_std(meanspec_thirteenco)) & (meanspec_thirteenco.spectral_axis>-200*u.km/u.s) & (meanspec_thirteenco.spectral_axis < 200 *u.km/u.s)
        thirteencomed = np.nanmedian(meanspec_thirteenco[~thirteencomask])
        meanspec_thirteenco_bl = meanspec_thirteenco - thirteencomed
        thirteenco_mom0 = u.Quantity(meanspec_thirteenco_bl[thirteencomask].sum() * dx_thirteenco)
        thirteenco_mom1 = u.Quantity((meanspec_thirteenco_bl[thirteencomask] * meanspec_thirteenco.spectral_axis[thirteencomask]).sum() / meanspec_thirteenco_bl[thirteencomask].sum())
        thirteenco_mom2 = u.Quantity((meanspec_thirteenco_bl[thirteencomask] * (meanspec_thirteenco.spectral_axis[thirteencomask]-thirteenco_mom1)**2).sum() / meanspec_thirteenco_bl[thirteencomask].sum())**0.5 * SIGMA2FWHM

        #If we want to fit the brightest line, we can use this
        sp_thirteenco = pyspeckit.Spectrum.from_hdu(meanspec_thirteenco.hdu)
        sp_thirteenco.xarr.convert_to_unit(u.km/u.s)
        #sp.plotter(figure=pl.figure(1))
        sp_thirteenco.specfit(guesses=[meanspec_thirteenco[thirteencoargmax].value, vmax.to(u.km/u.s).value, thirteenco_mom2.to(u.km/u.s).value])

        new_columns['13CO_mom0'].append(thirteenco_mom0.to(u.K*u.km/u.s).value)
        #new_columns['13CO_mom1'].append(thirteenco_mom1.to(u.km/u.s).value)
        #new_columns['13CO_fwhm'].append(thirteenco_mom2.to(u.km/u.s).value)
        new_columns['13CO_mom1'].append(sp_thirteenco.specfit.parinfo.SHIFT0.value)
        new_columns['13CO_fwhm'].append(sp_thirteenco.specfit.parinfo.WIDTH0.value * SIGMA2FWHM)
        new_columns['13CO_spectrum'].append(meanspec_thirteenco)
    else:
        new_columns['13CO_mom0'].append(np.nan)
        new_columns['13CO_mom1'].append(np.nan)
        new_columns['13CO_fwhm'].append(np.nan)
        new_columns['13CO_spectrum'].append(None)


    new_columns['has_zone_outside_data'].append(has_zone_outside_data)
    new_columns['npix_total'].append(npix)
    npix_cubes = submask.sum()
    new_columns['npix_cubes'].append(npix_cubes)

    if structure.idx in (8,9,11):
        assert npix_cubes != npix

    #if structure.idx == 50:
    #    print("HCN: {0},{1},{2}".format(hcn_mom0,hcn_mom1,hcn_mom2))
    #    print("HC3N: {0},{1},{2}".format(hc3n_mom0,hc3n_mom1,hc3n_mom2))
    #    import ipdb
    #    ipdb.set_trace()



cat = Table.read('column_and_temperature_catalog.ipac', format='ascii.ipac')

for key in units:
    cat.add_column(Column(data=new_columns[key], name=key, unit=units[key]))

cat.write('column_and_temperature_withHCN_catalog.ipac', format='ascii.ipac', overwrite=True)

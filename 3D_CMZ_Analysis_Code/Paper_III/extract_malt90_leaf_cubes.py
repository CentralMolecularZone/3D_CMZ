"""
Loops over HNCO, HCN, & HC3N maps from the MALT90 survey to extract sub-cubes.
Each sub-cube corresponds to a leaf (cloud) from our dendrogram.
Also inverts each masked sub-cube. These inverted cubes are used to obtain
averaged spectra from the medium around each cloud, which will be used to
perform a rough background subtraction for the spectra.
"""
from __future__ import division
import os
import numpy as np
import astrodendro
from astrodendro import Dendrogram, pp_catalog
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from spectral_cube import SpectralCube
import spectral_cube
# _____________________________________
file = "./../Data/higal_data/column_properunits_conv36_source_only.fits"
hdu = fits.open(file)[0]

dend = astrodendro.Dendrogram.load_from("./../Dendrogram_files/clouds_only_dendrogram.fits")
leaves = dend.leaves[9:(len(dend.leaves)-3)]

colfile = file
header = fits.getheader(colfile)
mywcs = wcs.WCS(header)

molecules = ["HNCO", "HCN", "HC3N"]

for mol in molecules:
    if mol == "HNCO":
        cube_file = './../Data/MALT90_data/CMZ_3mm_HNCO.fits'
    if mol == "HCN":
        cube_file = './../Data/MALT90_data/CMZ_3mm_HCN.fits'
    if mol == "HC3N":
        cube_file = './../Data/MALT90_data/CMZ_3mm_HC3N.fits'

    cube = SpectralCube.read(cube_file)
    cube_header = cube.header.copy()
    cube_header.update(mywcs.to_header())
    cube_header['NAXIS1'] = header['NAXIS1']
    cube_header['NAXIS2'] = header['NAXIS2']

    reproj_cube = cube.reproject(cube_header)

    i = 1
    for structure_indx in range(len(leaves)):
        if structure_indx == 45:
            continue
        else:
            structure = leaves[structure_indx]
            leaf_obj_mask = structure.get_mask()
            npix = leaf_obj_mask.sum()
            leaf_inds = structure.indices()
            view = [slice(leaf_inds[0].min(), leaf_inds[0].max()+1),
                    slice(leaf_inds[1].min(), leaf_inds[1].max()+1)]
            submask = leaf_obj_mask[view]
            submask_inv = ~submask
            cubeview = [slice(None),
                    slice(leaf_inds[0].min(), leaf_inds[0].max()+1),
                    slice(leaf_inds[1].min(), leaf_inds[1].max()+1)]

            cropcube = reproj_cube[cubeview].with_mask(submask[None,:,:])
            cropcube.write("./Leaf_cubes_MALT90/"+str(i)+"_"+str(mol)+"_cube.fits",
                           overwrite=True)

            cropcube_inv = reproj_cube[cubeview].with_mask(submask_inv[None,:,:])
            cropcube_inv.write("./Leaf_cubes_MALT90/Inverted/"+str(i)+"_"+str(mol)+"_cube.fits",
                               overwrite=True)
            i=i+1

    dend_brick = astrodendro.Dendrogram.load_from("./../Dendrogram_files/separate_brick16_dendrogram.fits")
    leaves_brick = dend_brick.leaves[-5]

    structure = leaves_brick
    leaf_obj_mask = structure.get_mask()
    npix = leaf_obj_mask.sum()
    leaf_inds = structure.indices()
    view = [slice(leaf_inds[0].min(), leaf_inds[0].max()+1),
            slice(leaf_inds[1].min(), leaf_inds[1].max()+1)]
    submask = leaf_obj_mask[view]
    submask_inv = ~submask
    cubeview = [slice(None),
            slice(leaf_inds[0].min(), leaf_inds[0].max()+1),
            slice(leaf_inds[1].min(), leaf_inds[1].max()+1)]

    cropcube = reproj_cube[cubeview].with_mask(submask[None,:,:])
    cropcube.write("./../Leaf_cubes_MALT90/15_"+str(mol)+"_cube.fits",
                   overwrite=True)

    cropcube_inv = reproj_cube[cubeview].with_mask(submask_inv[None,:,:])
    cropcube_inv.write("./../Leaf_cubes_MALT90/Inverted/15_"+str(mol)+"_cube.fits",
                       overwrite=True)

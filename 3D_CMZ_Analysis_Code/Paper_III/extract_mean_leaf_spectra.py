"""
Extracts spectra averaged over each leaf for all MALT90 & APEX sub-cubes.
"""
from __future__ import division
import os
import glob
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from spectral_cube import SpectralCube
import spectral_cube
# _____________________________________

# MALT90
os.chdir("./../Leaf_cubes_MALT90/")
if os.path.isdir("./meanspec/") == False:
    os.mkdir("./meanspec/")

for filename in glob.glob("*.fits"):
    t = os.path.splitext(filename)
    f = str(t[0])
    cube = SpectralCube.read(filename)
    meanspec = cube.mean(axis=(1,2))
    assert meanspec.size == cube.shape[0]
    meanspec.write("./meanspec/"+f+"_meanspec.fits", overwrite=True)

os.chdir("./../Leaf_cubes_MALT90/Inverted")
if os.path.isdir("./meanspec/") == False:
    os.mkdir("./meanspec/")

for filename in glob.glob("*.fits"):
    t = os.path.splitext(filename)
    f = str(t[0])
    cube = SpectralCube.read(filename)
    meanspec = cube.mean(axis=(1,2))
    assert meanspec.size == cube.shape[0]
    meanspec.write("./meanspec/"+f+"_meanspec.fits", overwrite=True)

# APEX
os.chdir("./../../Leaf_cubes_APEX/")
if os.path.isdir("./meanspec/") == False:
    os.mkdir("./meanspec/")

for filename in glob.glob("*.fits"):
    t = os.path.splitext(filename)
    f = str(t[0])
    cube = SpectralCube.read(filename)
    meanspec = cube.mean(axis=(1,2))
    assert meanspec.size == cube.shape[0]
    meanspec.write("./meanspec/"+f+"_meanspec.fits", overwrite=True)

os.chdir("./../Leaf_cubes_APEX/Inverted/")
if os.path.isdir("./meanspec/") == False:
    os.mkdir("./meanspec/")

for filename in glob.glob("*.fits"):
    t = os.path.splitext(filename)
    f = str(t[0])
    cube = SpectralCube.read(filename)
    meanspec = cube.mean(axis=(1,2))
    assert meanspec.size == cube.shape[0]
    meanspec.write("./meanspec/"+f+"_meanspec.fits", overwrite=True)

os.chdir("./../../Scripts")

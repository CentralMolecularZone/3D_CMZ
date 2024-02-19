"""
Plots a 6-panel figure displaying spectra and background-subtracted spectra
for each dendrogram leaf and all molecular line tracers in our sample.
"""
import os
import re
import numpy as np
import pylab as pl
import glob
import pyspeckit
from astropy.table import Table
from astropy import units as u
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)

for i in range(1,23):

    if i==3:
        continue
    else:
        fig = pl.figure(1)
        fig.clf()

        ax = pl.subplot(3,2,1)
        sp = pyspeckit.Spectrum("./../Leaf_cubes_MALT90/meanspec/"+str(i)+"_HCN_cube_meanspec.fits")
        sp_inv = pyspeckit.Spectrum("./../Leaf_cubes_MALT90/Inverted/meanspec/"+str(i)+"_HCN_cube_meanspec.fits")
        sp2 = sp - sp_inv
        sp2.plotter(figure=fig, axis=ax, clear=False, linestyle='--')
        sp.plotter(figure=fig, axis=ax, clear=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.set_ylabel("")
        ax.text(0.02, 0.95, "HCN",
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=8)
        ax.set_xlim(-200,200)

        ax = pl.subplot(3,2,2)
        sp = pyspeckit.Spectrum("./../Leaf_cubes_MALT90/meanspec/"+str(i)+"_HC3N_cube_meanspec.fits")
        sp_inv = pyspeckit.Spectrum("./../Leaf_cubes_MALT90/Inverted/meanspec/"+str(i)+"_HC3N_cube_meanspec.fits")
        sp2 = sp - sp_inv
        sp2.plotter(figure=fig, axis=ax, clear=False, linestyle='--')
        sp.plotter(figure=fig, axis=ax, clear=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.set_ylabel("")
        ax.text(0.02, 0.95, "HC$_{3}$N",
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=8)
        ax.annotate("-  Data", xy=(0.77, 0.92), xycoords='axes fraction', size=6)
        ax.annotate("-- Data - BG", xy=(0.77, 0.84), xycoords='axes fraction', size=6)
        ax.set_xlim(-200,200)

        ax = pl.subplot(3,2,3)
        sp = pyspeckit.Spectrum("./../Leaf_cubes_MALT90/meanspec/"+str(i)+"_HNCO_cube_meanspec.fits")
        sp_inv = pyspeckit.Spectrum("./../Leaf_cubes_MALT90/Inverted/meanspec/"+str(i)+"_HNCO_cube_meanspec.fits")
        sp2 = sp - sp_inv
        sp2.plotter(figure=fig, axis=ax, clear=False, linestyle='--')
        sp.plotter(figure=fig, axis=ax, clear=False)
        ax.set_ylabel("")
        ax.axes.get_xaxis().set_visible(False)
        ax.text(0.02, 0.95, "HNCO",
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=8)
        ax.set_xlim(-200,200)

        ax = pl.subplot(3,2,4)
        sp = pyspeckit.Spectrum("./../Leaf_cubes_APEX/meanspec/"+str(i)+"_13CO_cube_meanspec.fits")
        sp_inv = pyspeckit.Spectrum("./../Leaf_cubes_APEX/Inverted/meanspec/"+str(i)+"_13CO_cube_meanspec.fits")
        sp2 = sp - sp_inv
        sp2.plotter(figure=fig, axis=ax, clear=False, linestyle='--')
        sp.plotter(figure=fig, axis=ax, clear=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.set_ylabel("")
        ax.text(0.02, 0.95, "$^{13}$CO",
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=8)
        ax.set_xlim(-200,200)

        ax = pl.subplot(3,2,5)
        sp = pyspeckit.Spectrum("./../Leaf_cubes_APEX/meanspec/"+str(i)+"_C18O_cube_meanspec.fits")
        sp_inv = pyspeckit.Spectrum("./../Leaf_cubes_APEX/Inverted/meanspec/"+str(i)+"_C18O_cube_meanspec.fits")
        sp2 = sp - sp_inv
        sp2.plotter(figure=fig, axis=ax, clear=False, linestyle='--')
        sp.plotter(figure=fig, axis=ax, clear=False)
        ax.xaxis.label.set_visible(False)
        ax.set_ylabel("")
        ax.text(0.02, 0.95, "C$^{18}$O",
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=8)
        ax.set_xlim(-200,200)

        ax = pl.subplot(3,2,6)
        sp = pyspeckit.Spectrum("./../Leaf_cubes_APEX/meanspec/"+str(i)+"_H2CO_303_202_cube_meanspec.fits")
        sp_inv = pyspeckit.Spectrum("./../Leaf_cubes_APEX/Inverted/meanspec/"+str(i)+"_H2CO_303_202_cube_meanspec.fits")
        sp2 = sp - sp_inv
        sp2.plotter(figure=fig, axis=ax, clear=False, linestyle='--')
        sp.plotter(figure=fig, axis=ax, clear=False)
        ax.set_ylabel("")
        ax.xaxis.label.set_visible(False)
        ax.text(0.02, 0.95, "H$_{2}$CO$_{303-202}$",
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=8)
        ax.set_xlim(-200,200)

        fig.text(0.5, 0.02, 'Velocity (km/s)', ha='center', fontsize=10)
        fig.text(0.04, 0.5, 'Brightness Temperature (K)', va='center', rotation='vertical', fontsize=10)
        sp.plotter.savefig("./../Figs/Multispec/bg-sub/"+str(i)+"_multispec.pdf")
        fig.clf()

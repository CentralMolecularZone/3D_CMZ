import os
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.table import Table
import matplotlib.pyplot as pl
import matplotlib.colors as mc
from scipy.interpolate import NearestNDInterpolator
from numpy import linspace, array, logspace, sin, cos, pi, arange, sqrt, arctan2, arccos
from mpl_toolkits.mplot3d import Axes3D
from coords import *
plt.style.use('classic')
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.size"] = "16"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

def rotate(x,y,theta):
  xprime = x*cos(theta) - y*sin(theta)
  yprime = x*sin(theta) + y*cos(theta)
  return xprime, yprime

#############################
# create class that stores a ring
#############################

class Ring:

    def __init__(self,t,a,b,z,v0,theta,xyzsun,vxyzsun):
        self.t     = t
        self.a     = a
        self.b     = b
        self.z0    = z
        self.v0    = v0
        self.theta = theta
        self.x     = a*cos(t)
        self.y     = -b*sin(t)
        self.z     = self.z0*np.ones(t.size)
        self.R     = sqrt(self.x**2+self.y**2)
        self.phi   = -arctan2(self.y,self.x)
        self.ephix = -sin(self.phi) # unit vector parallel to circle
        self.ephiy = -cos(self.phi) # unit vector parallel to circle
        norm       = sqrt((a*sin(t))**2+(b*cos(t))**2)
        self.ex    = -a*sin(t)/norm # unit vector parellel to ellipse
        self.ey    = -b*cos(t)/norm  # unit vector parallel to ellipse
        self.cosalpha = self.ex*self.ephix + self.ey*self.ephiy
        self.vphi  = self.R[0]*self.v0/self.R # assume conservation of angular momentum
        self.v     = self.vphi/self.cosalpha # total speed along the orbit
        self.vx    = +self.v*self.ex
        self.vy    = +self.v*self.ey
        self.vz    = np.zeros(t.size)
        self.x,self.y   = rotate(self.x,self.y,theta)
        self.vx,self.vy = rotate(self.vx,self.vy,theta)
        self.X,self.Y,self.Z,self.vX,self.Vy,self.vZ = xyz2XYZ(self.x,self.y,self.z,self.vx,self.vy,self.vz,xyzsun[0],xyzsun[1],xyzsun[2],vxyzsun[0],vxyzsun[1],vxyzsun[2])
        self.l,self.b,self.r,self.vl,self.vb,self.vr = xyz2lbr(self.x,self.y,self.z,self.vx,self.vy,self.vz,xyzsun[0],xyzsun[1],xyzsun[2],vxyzsun[0],vxyzsun[1],vxyzsun[2])
        self.mu_l, self.mu_b = vlb_2_mulb(self.r,self.vl*100,self.vb*100)

#############################
# define sun position & velocity
#############################

xsun  = 0.0
ysun  = -8.0
zsun  = 0.0

vxsun = -2.2
vysun = 0.0
vzsun = 0.0

xyzsun  = [xsun, ysun, zsun ]
vxyzsun = [vxsun,vysun,vzsun]

phisun  = arctan2(ysun,xsun)

#############################
# define rings
#############################

a_rings = np.array([0.1,0.13,0.2])*1 # x axis of the rings in kpc
b_rings = np.array([0.1,0.13,0.2])*1.5 # y axis of the rings in kpc
z_rings = np.array([0.0,0.0,0.0] ) # height of the ring above/below plane in kpc
v_rings = np.array([1.0,1.0,1.0] ) # tangential velocity at initial point
theta_rings = [np.radians(20),np.radians(20),np.radians(20)] # inclination of the rings

dphi = np.radians(1)
phi  = arange(phisun+dphi/2,phisun+2*pi,dphi)

Rings = [Ring(phi,a_rings[i],b_rings[i],z_rings[i],v_rings[i]*100,theta_rings[i],xyzsun,vxyzsun) for i in arange(len(a_rings))]

scouse_fits = pd.read_csv('final_cmz_scouse_hnco_fits.csv',
              usecols=[0,1,2,3,5,7], names=['n', 'l', 'b', 'amp', 'velocity', 'FWHM'],
              sep="\s+")

cat = pd.read_csv('CMZ_cloud_catalogue_data.csv')

absorp_tab = Table.read('absorption_comparison_table.tex')
absorp_value = absorp_tab['fraction_value']
absorp_13co = absorp_tab['co_h2co_ratio']
absorp_index = absorp_tab['leaf_id']
divnorm = colors.TwoSlopeNorm(vmin=absorp_value.min(), vcenter=1.0, vmax=absorp_value.max())
divnorm_13co = colors.TwoSlopeNorm(vmin=absorp_13co.min(), vcenter=1.0, vmax=absorp_13co.max())

l   = scouse_fits['l'].values
b   = scouse_fits['b'].values
v   = scouse_fits['velocity'].values
amp = scouse_fits['amp'].values

fig, ax = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(box_aspect=1), figsize=(15,15))
plt.subplots_adjust(wspace=0, hspace=0)

for column in ax:
    for row in column:
        row.scatter(l, v, c=amp, cmap='Blues', marker='o', s=1, facecolor='1',
                    norm=mpl.colors.LogNorm(), alpha=0.2, zorder=1)
        for i in range(0,len(cat)):
            if cat['index'][i] in absorp_index:
                sc = row.scatter(cat['l'][i], cat['v'][i], marker='o', s=100, edgecolor='k', zorder=3, c=absorp_value[absorp_index == cat['index'][i]], cmap='bwr', norm=divnorm)
            else:
                row.scatter(cat['l'][i], cat['v'][i], marker='*', s=150, edgecolor='k', zorder=3, c='k')
        row.grid(True)
        row.set_xlim(1,-0.6)
        row.set_ylim(-140,140)
        row.set_rasterized(True)

"""
KDL
"""
kdl = np.loadtxt('kdl_lbv.txt', unpack=True)
lk, bk, vk = kdl[0], kdl[1], kdl[2]
ax[0][0].plot(lk[0:16], vk[0:16], ls='-', color='blue', linewidth=2, zorder=2)
ax[0][0].plot(lk[15:26], vk[15:26], ls='-', color='red', linewidth=2, zorder=2, alpha=0.4)
ax[0][0].plot(lk[25:36], vk[25:36], ls='-', color='red', linewidth=2, zorder=2, alpha=0.4)
ax[0][0].plot(lk[35:51], vk[35:51], ls='-', color='blue', linewidth=2, zorder=2)
ax[0][0].xaxis.set_ticklabels([])
ax[0][0].yaxis.set_ticklabels([])
ax[0][0].text(0.94, 0.96,'KDL', ha='center', va='center', transform=ax[0][0].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

"""
Sofue
"""
sofue1 = np.loadtxt('sofue_lbv_arm1.txt', unpack=True)
sofue1_l, sofue1_b, sofue1_v = sofue1[0], sofue1[1], sofue1[2]
ax[0][1].plot(sofue1_l, sofue1_v, ls='-', color='blue', linewidth=2, zorder=2)

sofue2 = np.loadtxt('sofue_lbv_arm2.txt', unpack=True)
sofue2_l, sofue2_b, sofue2_v = sofue2[0], sofue2[1], sofue2[2]
ax[0][1].plot(sofue2_l, sofue2_v, ls='-', color='red', linewidth=2, zorder=2, alpha=0.4)
ax[0][1].xaxis.set_ticklabels([])
ax[0][1].yaxis.set_ticklabels([])
ax[0][1].text(0.92, 0.96,'Sofue', ha='center', va='center', transform=ax[0][1].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

"""
Molinari
"""
ml = np.loadtxt('molinari_lbvlos.txt', unpack=True)
lk, bk, vk = ml[0], ml[1], ml[2]
ax[1][0].plot(-1*lk[80:580], vk[80:580]-14, ls='-', color='red', linewidth=2, zorder=2, alpha=0.4)
ax[1][0].plot(-1*lk[0:80], vk[0:80]-14, ls='-', color='blue', linewidth=2, zorder=2)
ax[1][0].plot(-1*lk[581:], vk[581:]-14, ls='-', color='blue', linewidth=2, zorder=2)
ax[1][0].text(0.89, 0.96,'Molinari', ha='center', va='center', transform=ax[1][0].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))
ax[1][0].set_xlabel('Galactic longitude')
ax[1][0].set_ylabel('Velocity (km/s)')

"""
Rings
"""
i = 0
colours = ['Green', 'Orange', 'Magenta']
for Ring in Rings:
    if Ring == Rings[1]:
        ax[1][1].plot(np.degrees(Ring.l),Ring.vr,c='black', linewidth=2)

ax[1][1].xaxis.set_ticklabels([])
ax[1][1].yaxis.set_ticklabels([])
ax[1][1].text(0.91, 0.96,'Ellipse', ha='center', va='center', transform=ax[1][1].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

cax = fig.add_axes([0.128, .92, 0.77, 0.03])
cb = plt.colorbar(sc, ax=ax.ravel().tolist(), orientation='horizontal', cax=cax)
cb.ax.set_title('Median 4.8 GHz continuum / Minimum H$_{2}$CO (1$_{1,0}$-1$_{1,1}$)', size=20, fontweight='bold', fontname='sans-serif')
plt.savefig('ellipses_lv_colorbar_H2CO.pdf')
plt.clf()

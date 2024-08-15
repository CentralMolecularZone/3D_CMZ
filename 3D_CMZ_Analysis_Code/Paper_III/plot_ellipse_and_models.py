from coords import *
import numpy as np
import pandas as pd
from astropy.table import Table
import matplotlib.pyplot as plt
from numpy import sin, cos, pi, arange, sqrt, arctan2
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.table import Table
import matplotlib.colors as mc
import matplotlib.patches as patches
from numpy import sin, cos, pi, arange, sqrt, arctan2
import matplotlib.patches as patches
import matplotlib.patheffects as pe
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

Rsun = 8.1
l_offset = np.radians(0.05)
b_offset = np.radians(-0.0462)
x_offset = Rsun*l_offset
y_offset = 0

class Ring:

    def __init__(self,t,a,b,z,v0,theta,xyzsun,vxyzsun,alpha=0.4):
        self.t     = t
        self.a     = a
        self.b     = b
        self.z0    = z
        self.v0    = v0
        self.theta = theta
        self.x     = a*cos(t)
        self.y     = -b*sin(t)
        self.z     = self.z0*sin(-2*t + alpha)
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
        l_offset = np.radians(0.05)
        b_offset = np.radians(-0.0462)
        x_offset = Rsun*l_offset
        y_offset = 0
        z_offset = Rsun*b_offset
        self.l,self.b,self.r,self.vl,self.vb,self.vr = xyz2lbr(
            self.x,self.y,self.z,self.vx,self.vy,self.vz,
            xyzsun[0],xyzsun[1],xyzsun[2],vxyzsun[0],vxyzsun[1],vxyzsun[2])
        self.l += l_offset
        self.b += b_offset
        self.x += x_offset
        self.y += y_offset
        self.z += z_offset
        self.mu_l, self.mu_b = vlb_2_mulb(self.r, self.vl*100, self.vb*100)
        self.mu_l, self.mu_b = vlb_2_mulb(self.r,self.vl*100,self.vb*100)

#############################

#############################
# define sun position & velocity
#############################

xsun  = 0.0
ysun  = -8.1
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

a_rings = 0.09 # x axis of the rings in kpc
b_rings = 0.055 # y axis of the rings in kpc
z_rings = 0.0125 # height of the ring above/below plane in kpc
v_rings = 130.0 # tangential velocity at initial point
theta_rings = np.radians(25) # rotation of the rings
dphi = np.radians(1)
phi  = arange(phisun+dphi/2,phisun+2*pi,dphi)
Rings = Ring(phi,a_rings,b_rings,z_rings,v_rings,theta_rings,xyzsun,vxyzsun)

def plot_ellipse_model(ring, x_offset, y_offset, l_offset, b_offset):
    """
    Plot the ellipse model.

    Args:
        ring (Ring): The ring object to plot.
        x_offset (float): The x offset.
        y_offset (float): The y offset.
        l_offset (float): The l offset.
        b_offset (float): The b offset.
    """
    fig, ax = plt.subplots(nrows=1, ncols=3, subplot_kw=dict(box_aspect=1), figsize=(20,6))
    plt.subplots_adjust(wspace=0.4, hspace=0)

    for i in range(len(ring.y) - 1):
        colour = 'red' if ring.y[i] > 0 else 'blue'
        linstyle = 'dotted' if ring.y[i] > 0 else 'dashed'

        ax[0].plot(ring.x[i:i+2], ring.y[i:i+2], color=colour, linewidth=10, alpha=0.05)
        ax[0].plot(ring.x[i:i+2], ring.y[i:i+2], color=colour, linewidth=1, ls=linstyle)
        ax[1].plot(np.degrees(ring.l[i:i+2]), np.degrees(ring.b[i:i+2]), c=colour, linewidth=10, alpha=0.05)
        ax[1].plot(np.degrees(ring.l[i:i+2]), np.degrees(ring.b[i:i+2]), c=colour, linewidth=1, ls=linstyle)
        ax[2].plot(np.degrees(ring.l[i:i+2]), ring.vr[i:i+2], c=colour, linewidth=10, alpha=0.05)
        ax[2].plot(np.degrees(ring.l[i:i+2]), ring.vr[i:i+2], c=colour, linewidth=1, ls=linstyle)

    ax[0].plot(x_offset, y_offset, marker='v', markersize=12, color='k')
    ax[0].plot(-x_offset, y_offset, marker='*', markersize=15, color='k')
    ax[0].plot([-x_offset, x_offset], [y_offset, y_offset], linestyle='dotted', color='k', linewidth=2)
    ax[0].grid(True)
    ax[0].set_ylabel('y (kpc)')
    ax[0].set_xlabel('x (kpc)')
    ax[0].set_xlim(-0.09, 0.10)
    ax[0].set_ylim(-0.075, 0.075)

    ax[1].plot(np.degrees(l_offset), np.degrees(b_offset), marker='v', markersize=12, color='black')
    ax[1].plot(-np.degrees(l_offset), np.degrees(b_offset), marker='*', markersize=15, color='black')
    ax[1].plot([-np.degrees(l_offset), np.degrees(l_offset)], [np.degrees(b_offset), np.degrees(b_offset)], linestyle='dotted', color='black', linewidth=2)
    ax[1].grid(True)
    ax[1].set_ylabel('$b$ (deg)')
    ax[1].set_xlabel('$l$ (deg)')
    ax[1].set_xlim(0.75, -0.65)
    ax[1].set_ylim(-0.15, 0.06)

    ax[2].plot(np.degrees(l_offset), 0, marker='v', markersize=12, color='k')
    ax[2].plot(-np.degrees(l_offset), 0, marker='*', markersize=15, color='k')
    ax[2].plot([-np.degrees(l_offset), np.degrees(l_offset)], [0, 0], linestyle='dotted', color='k', linewidth=2)
    ax[2].grid(True)
    ax[2].set_ylabel('$v$ (km s$^{-1}$)')
    ax[2].set_xlabel('$l$ (deg)')
    ax[2].set_xlim(0.75, -0.65)
    ax[2].set_ylim(-120, 120)

    plt.savefig('ellipse_model.pdf')
    plt.clf()
    plt.close()

plot_ellipse_model(Rings, x_offset, y_offset, l_offset, b_offset)

scouse_fits = pd.read_csv('final_cmz_scouse_hnco_fits.csv',
              usecols=[0,1,2,3,5,7], names=['n', 'l', 'b', 'amp', 'velocity', 'FWHM'],
              sep=r"\s+")

cat = pd.read_csv('full_catalogue.csv')

absorp_tab = Table.read('absorption_comparison_table_velocity_restricted_weighted.tex')
absorp_value = absorp_tab['fraction_value']
absorp_index = absorp_tab['leaf_id']
divnorm = mc.TwoSlopeNorm(vmin=absorp_value.min(), vcenter=1.0, vmax=2.0)

l   = scouse_fits['l'].values
b   = scouse_fits['b'].values
v   = scouse_fits['velocity'].values
amp = scouse_fits['amp'].values

fig, ax = plt.subplots(nrows=2, ncols=2,  figsize=(20,15))
plt.subplots_adjust(wspace=0, hspace=0)

base_marker_size = 100
normalized_masses = cat['mass'].values / 2e4
marker_sizes = normalized_masses * base_marker_size
    
for column in ax:
    for row in column:
        for i in range(0,len(cat)):
            size = marker_sizes[i]
            ref_marker_position = (0.05, 0.05)
            ref_marker_size = 100
            ax[1][1].scatter(0.05, 0.95, s=ref_marker_size, color='white', edgecolor='black', linewidth=2, 
                            transform=ax[1][1].transAxes, zorder=4)
            ax[1][1].annotate(" = 2$\\times$10$^{4}$ M$_{\odot}$", 
                            xy=(0.05, 0.95), xycoords='axes fraction',
                            xytext=(5, 0), textcoords='offset points',
                            ha='left', va='center', fontsize=16, color='black', zorder=5, fontfamily='sans-serif')

            indices = np.where(absorp_index == cat['index'][i])[0]
            if cat['index'][i] in absorp_index:
                for idx in indices:
                    if absorp_tab['abs'][idx] == 0:
                        colour = 'red'
                        hat = ''
                    elif absorp_tab['abs'][idx] == 1:
                        colour = 'blue'
                        hat = ''
                    elif absorp_tab['abs'][idx] == 2 and not cat['index'][i] == 8:
                        colour = 'lightgrey'
                        hat = '//'
                    row.scatter(cat['l'][i], cat['b'][i], marker='o', s=size, edgecolor='k', zorder=3, c=colour, hatch=hat)
            else:
                row.scatter(cat['l'][i], cat['b'][i], marker='s', s=size, edgecolor='k', zorder=3, c='None')

            marker_coords = (cat['l'][i], cat['b'][i])
            source_index = cat['index'][i]
            if not source_index >=29 and not source_index == 21 and not source_index == 14:
                ax[1][1].text(marker_coords[0], marker_coords[1]-0.01, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
            if source_index == 21:
                ax[1][1].text(marker_coords[0], marker_coords[1]+0.01, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
            if source_index == 14:
                ax[1][1].text(marker_coords[0], marker_coords[1]+0.01, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
        ax[1][1].plot(np.degrees(l_offset), np.degrees(b_offset), marker='v', markersize=12, color='black')
        row.plot(-np.degrees(l_offset), np.degrees(b_offset), marker='*', markersize=15, color='black')
        ax[1][1].plot([-np.degrees(l_offset), np.degrees(l_offset)], [np.degrees(b_offset), np.degrees(b_offset)], linestyle='dotted', color='black', linewidth=2)
        row.grid(True)
        row.set_xlim(1,-0.65)
        row.set_ylim(-0.25,0.1)
        row.set_rasterized(True)


x_start = 0.7
x_end = 0.8
y_start = 0
y_end = -0.1
rect = patches.Rectangle((x_start, y_start), x_end - x_start, y_end - y_start, linewidth=1, linestyle='dashed',
                         edgecolor='red', facecolor='none')
"""
KDL
"""
kdl = np.loadtxt('kdl_lbv.txt', unpack=True)
lk, bk, vk = kdl[0], kdl[1], kdl[2]
ax[0][0].plot(lk[0:16], bk[0:16], ls='-', color='blue', alpha=0.2, linewidth=10)
ax[0][0].plot(lk[0:16], bk[0:16], ls='-', color='blue', linewidth=1)
ax[0][0].plot(lk[15:26], bk[15:26], ls='-', color='red', alpha=0.2, linewidth=10)
ax[0][0].plot(lk[15:26], bk[15:26], ls='dotted', color='red', linewidth=1)
ax[0][0].plot(lk[25:36], bk[25:36], ls='-', color='red', alpha=0.2, linewidth=10)
ax[0][0].plot(lk[25:36], bk[25:36], ls='dotted', color='red', linewidth=1)
ax[0][0].plot(lk[35:51], bk[35:51], ls='-', color='blue', alpha=0.2, linewidth=10)
ax[0][0].plot(lk[35:51], bk[35:51], ls='-', color='blue', linewidth=1)
ax[0][0].xaxis.set_ticklabels([])
ax[0][0].yaxis.set_ticklabels([])
ax[0][0].text(0.943, 0.96,'KDL', ha='center', va='center', transform=ax[0][0].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

"""
Sofue
"""
sofue1 = np.loadtxt('sofue_lbv_arm1.txt', unpack=True)
sofue1_l, sofue1_b, sofue1_v = sofue1[0], sofue1[1], sofue1[2]
ax[0][1].plot(sofue1_l, sofue1_b, ls='-', color='blue', alpha=0.2, linewidth=10)
ax[0][1].plot(sofue1_l, sofue1_b, ls='-', color='blue', linewidth=1)

sofue2 = np.loadtxt('sofue_lbv_arm2.txt', unpack=True)
sofue2_l, sofue2_b, sofue2_v = sofue2[0], sofue2[1], sofue2[2]
ax[0][1].plot(sofue2_l, sofue2_b, ls='-', color='red', linewidth=10, alpha=0.2)
ax[0][1].plot(sofue2_l, sofue2_b, ls='dotted', color='red', linewidth=1)
ax[0][1].xaxis.set_ticklabels([])
ax[0][1].yaxis.set_ticklabels([])
ax[0][1].text(0.93, 0.96,'Sofue', ha='center', va='center', transform=ax[0][1].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

"""
Molinari
"""
ml = np.loadtxt('molinari_lbvlos.txt', unpack=True)
lk, bk, vk = ml[0], ml[1], ml[2]
ax[1][0].plot(-1*lk[80:580], bk[80:580], ls='-', color='red', linewidth=10, alpha=0.2)
ax[1][0].plot(-1*lk[80:580], bk[80:580], ls='dotted', color='red', linewidth=1)
ax[1][0].plot(-1*lk[0:80], bk[0:80], ls='-', color='blue', alpha=0.2, linewidth=10)
ax[1][0].plot(-1*lk[0:80], bk[0:80], ls='-', color='blue', linewidth=1)
ax[1][0].plot(-1*lk[581:], bk[581:], ls='-', color='blue', alpha=0.2, linewidth=10)
ax[1][0].plot(-1*lk[581:], bk[581:], ls='-', color='blue', linewidth=1)
ax[1][0].text(0.91, 0.96,'Molinari', ha='center', va='center', transform=ax[1][0].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))
ax[1][0].set_xlabel('Galactic longitude')
ax[1][0].set_ylabel('Galactic latitude')

"""
Rings
"""
for i in range(len(Rings.y) - 1):
    colour = 'red' if Rings.y[i] > 0 else 'blue'
    linstyle = 'dotted' if Rings.y[i] > 0 else 'dashed'

    ax[1][1].plot(np.degrees(Rings.l[i:i+2]),np.degrees(Rings.b[i:i+2]),c=colour, linewidth=10, alpha=0.05)
    ax[1][1].plot(np.degrees(Rings.l[i:i+2]),np.degrees(Rings.b[i:i+2]),c=colour, linewidth=1, ls=linstyle)

ax[1][1].xaxis.set_ticklabels([])
ax[1][1].yaxis.set_ticklabels([])
ax[1][1].text(0.92, 0.96,'Ellipse', ha='center', va='center', transform=ax[1][1].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

plt.savefig('3d_models_lb.pdf')
plt.clf()
plt.close()

#############################
# Plot lv
#############################

cat = pd.read_csv('CMZ_cloud_catalogue_data.csv')
sigma = cat['sigma']
radius = cat['rad']
radius_deg = np.degrees(radius / 8100)

absorp_tab = Table.read('absorption_comparison_table_velocity_restricted_weighted.tex')
absorp_value = absorp_tab['fraction_value']
absorp_index = absorp_tab['leaf_id']
absorp = absorp_tab['abs']
vel          = absorp_tab['velocity']

fig, ax = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(box_aspect=1), figsize=(15,15))
plt.subplots_adjust(wspace=0, hspace=0)

for column in ax:
    for row in column:
        row.scatter(l, v, c=amp, cmap='Blues', marker='o', s=1, facecolor='1',
                    norm=mpl.colors.LogNorm(), alpha=0.2, zorder=1)
        for i in range(0,len(cat)):
            indices = np.where(absorp_index == cat['index'][i])[0]
            if cat['index'][i] in absorp_index:
                for idx in indices:
                    if absorp[(vel == cat['v'][i]) & (absorp_index == cat['index'][i])] == 0:
                        colour = 'red'
                        hat = ''
                    elif absorp[(vel == cat['v'][i]) & (absorp_index == cat['index'][i])] == 1:
                        colour = 'blue'
                        hat = ''
                    elif absorp[(vel == cat['v'][i]) & (absorp_index == cat['index'][i])] == 2:
                        colour = 'lightgrey'
                        hat = '//'
                    ellipse = patches.Ellipse((cat['l'][i], cat['v'][i]), width=radius_deg[i]*2, height=sigma[i], edgecolor='k', facecolor=colour, zorder=3, hatch=hat)
                    row.add_patch(ellipse)
            else:
                row.scatter(cat['l'][i], cat['v'][i], marker='s', s=150, edgecolor='k', zorder=1, c='white', alpha=0.65)

            marker_coords = (cat['l'][i], cat['v'][i])
            source_index = cat['index'][i]
            if not source_index >=29 and not source_index == 10 and not source_index == 9 and not source_index == 13 and not source_index == 16 and not source_index == 8 and not source_index == 6:
                ax[1][1].text(marker_coords[0], marker_coords[1]-8, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
            if source_index == 10:
                ax[1][1].text(marker_coords[0]-0.05, marker_coords[1]-8, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
            if source_index == 9:
                ax[1][1].text(marker_coords[0]+0.025, marker_coords[1]-8, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
            if source_index == 13:
                ax[1][1].text(marker_coords[0], marker_coords[1]+8, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
            if source_index == 16:
                ax[1][1].text(marker_coords[0], marker_coords[1]+8, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
            if source_index == 8:
                ax[1][1].text(marker_coords[0], marker_coords[1]+8, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])
            if source_index == 6:
                ax[1][1].text(marker_coords[0]-0.03, marker_coords[1]-2, str(source_index), 
                        ha='center', va='center', fontsize=12, color='white', path_effects=[pe.withStroke(linewidth=2,foreground="black")])


        ax[1][1].plot(np.degrees(l_offset), 0, marker='v', markersize=12, color='black')
        row.plot(-np.degrees(l_offset), 0, marker='*', markersize=15, color='black')
        ax[1][1].plot([-np.degrees(l_offset), np.degrees(l_offset)], [0,0], linestyle='dotted', color='black', linewidth=2)
        row.grid(True)
        row.set_xlim(1,-0.6)
        row.set_ylim(-140,140)
        row.set_rasterized(True)

"""
KDL
"""
kdl = np.loadtxt('kdl_lbv.txt', unpack=True)
lk, bk, vk = kdl[0], kdl[1], kdl[2]
ax[0][0].plot(lk[0:16], vk[0:16], ls='-', color='blue', alpha=0.2, linewidth=10)
ax[0][0].plot(lk[0:16], vk[0:16], ls='-', color='blue', linewidth=1)
ax[0][0].plot(lk[15:26], vk[15:26], ls='-', color='red', linewidth=10, alpha=0.2)
ax[0][0].plot(lk[15:26], vk[15:26], ls='dotted', color='red', linewidth=1)
ax[0][0].plot(lk[25:36], vk[25:36], ls='-', color='red', linewidth=10, alpha=0.2)
ax[0][0].plot(lk[25:36], vk[25:36], ls='dotted', color='red', linewidth=1)
ax[0][0].plot(lk[35:51], vk[35:51], ls='-', color='blue', alpha=0.2, linewidth=10)
ax[0][0].plot(lk[35:51], vk[35:51], ls='-', color='blue', linewidth=1)
ax[0][0].xaxis.set_ticklabels([])
ax[0][0].yaxis.set_ticklabels([])
ax[0][0].text(0.94, 0.96,'KDL', ha='center', va='center', transform=ax[0][0].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

"""
Sofue
"""
sofue1 = np.loadtxt('sofue_lbv_arm1.txt', unpack=True)
sofue1_l, sofue1_b, sofue1_v = sofue1[0], sofue1[1], sofue1[2]
ax[0][1].plot(sofue1_l, sofue1_v, ls='-', color='blue', alpha=0.2, linewidth=10)
ax[0][1].plot(sofue1_l, sofue1_v, ls='-', color='blue', linewidth=1)

sofue2 = np.loadtxt('sofue_lbv_arm2.txt', unpack=True)
sofue2_l, sofue2_b, sofue2_v = sofue2[0], sofue2[1], sofue2[2]
ax[0][1].plot(sofue2_l, sofue2_v, ls='-', color='red', linewidth=10, alpha=0.2)
ax[0][1].plot(sofue2_l, sofue2_v, ls='dotted', color='red', linewidth=1)
ax[0][1].xaxis.set_ticklabels([])
ax[0][1].yaxis.set_ticklabels([])
ax[0][1].text(0.92, 0.96,'Sofue', ha='center', va='center', transform=ax[0][1].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

"""
Molinari
"""
ml = np.loadtxt('molinari_lbvlos.txt', unpack=True)
lk, bk, vk = ml[0], ml[1], ml[2]
ax[1][0].plot(-1*lk[80:580], vk[80:580]-14, ls='-', color='red', linewidth=10, alpha=0.2)
ax[1][0].plot(-1*lk[80:580], vk[80:580]-14, ls='dotted', color='red', linewidth=1)
ax[1][0].plot(-1*lk[0:80], vk[0:80]-14, ls='-', color='blue', alpha=0.2, linewidth=10)
ax[1][0].plot(-1*lk[0:80], vk[0:80]-14, ls='-', color='blue', linewidth=1)
ax[1][0].plot(-1*lk[581:], vk[581:]-14, ls='-', color='blue', alpha=0.2, linewidth=10)
ax[1][0].plot(-1*lk[581:], vk[581:]-14, ls='-', color='blue', linewidth=1)
ax[1][0].text(0.89, 0.96,'Molinari', ha='center', va='center', transform=ax[1][0].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))
ax[1][0].set_xlabel('Galactic longitude')
ax[1][0].set_ylabel('Velocity (km/s)')

"""
Rings
"""
for i in range(len(Rings.y) - 1):
    colour = 'red' if Rings.y[i] > 0 else 'blue'
    linstyle = 'dotted' if Rings.y[i] > 0 else 'dashed'

    ax[1][1].plot(np.degrees(Rings.l[i:i+2]),Rings.vr[i:i+2],c=colour, linewidth=10, alpha=0.05)
    ax[1][1].plot(np.degrees(Rings.l[i:i+2]),Rings.vr[i:i+2],c=colour, linewidth=1, ls=linstyle)

ax[1][1].xaxis.set_ticklabels([])
ax[1][1].yaxis.set_ticklabels([])
ax[1][1].text(0.91, 0.96,'Ellipse', ha='center', va='center', transform=ax[1][1].transAxes,
           bbox=dict(facecolor='white', edgecolor='k', alpha=1))

plt.savefig('3d_models_lv.pdf')
plt.clf()
plt.close()
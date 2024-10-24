from __future__ import division
import numpy as np
import pandas as pd
from astropy import units as u
from astropy import constants as const
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib as mpl

def load_and_clean_data(kinematic_file, physical_file, drop_leaves):
    """Load and clean kinematic and physical data for CMZ clouds."""
    # Load kinematic data
    kin_data = Table.read(kinematic_file).to_pandas()
    kin_data = kin_data[~kin_data['Leaf ID'].isin(drop_leaves)]
    kin_data = kin_data[['Leaf ID', 'HNCO mom1', 'HNCO fwhm']].reset_index(drop=True).drop(labels=0)
    
    # Remove empty rows
    invalid_idx = kin_data.index[kin_data['HNCO mom1'] == '-'].tolist()
    kin_data = kin_data.drop(labels=invalid_idx).reset_index(drop=True)
    
    # Load physical data
    phys_cols = ['Leaf ID', 'Area', 'l', 'b', 'Median column density', 
                 'Mass', 'Radius', 'Median dust temperature']
    phys_data = Table.read(physical_file).to_pandas()
    phys_data = phys_data[~phys_data['Leaf ID'].isin(drop_leaves)]
    phys_data = phys_data[phys_cols].reset_index(drop=True).drop(labels=0)
    phys_data = phys_data.drop(labels=invalid_idx).reset_index(drop=True)
    
    return kin_data, phys_data

def calculate_properties(phys_data, kin_data):
    """Calculate surface density and sigma^2/R for CMZ clouds."""
    mass_values = phys_data['Mass'].values.astype(float)
    radius_values = phys_data['Radius'].values.astype(float)
    area_values = phys_data['Area'].values.astype(float)
    fwhm_values = kin_data['HNCO fwhm'].values.astype(float)
    
    Mass = mass_values * u.Msun
    Radius = radius_values * u.pc
    Area = area_values * (u.pc * u.pc)
    sigma = np.around((fwhm_values * (u.km/u.s)) / (2 * np.sqrt(2*np.log(2))), 1)
    
    surf_den = Mass.to(u.g) / Area.to(u.cm*u.cm)
    sigma2_r = np.around((sigma**2 / Radius), 1)
    
    return np.log10(surf_den.value), np.log10(sigma2_r.value)

def load_comparison_data(grs_file, f18_file, rosolowsky_file):
    """Load comparison data from other surveys."""
    # GRS data
    grs = np.loadtxt(grs_file)
    grs_x, grs_y = grs[:,0], grs[:,1]
    
    # Faesi et al. (2018) data
    f18 = np.loadtxt(f18_file)
    f18_mass_per_area = f18[:,0] * (u.Msun/(u.pc*u.pc))
    f18_x = np.log10(f18_mass_per_area.to(u.g/(u.cm*u.cm)).value)
    f18_y = np.log10(f18[:,1])
    
    # Rosolowsky et al. (2021) data
    r21 = Table.read(rosolowsky_file).to_pandas()
    r21 = r21[['SIGV_KMS', 'RAD_PC', 'MLUM_MSUN']]
    # 3D radius conversion criterion from Rosolowsky et al. (2021)
    r21['RAD_PC'] = np.where(r21['RAD_PC'] > 50, 
                            np.cbrt((r21['RAD_PC']**2 * 100) / 2), 
                            r21['RAD_PC'])
    
    r21_sigma = r21['SIGV_KMS'].values**2 / r21['RAD_PC'].values
    mass_msun = r21['MLUM_MSUN'].values * u.Msun
    radius_pc = r21['RAD_PC'].values * u.pc
    area_cm2 = (2*np.pi * radius_pc.to(u.cm)**2)
    r21_surf_den = mass_msun.to(u.g) / area_cm2
    
    return (grs_x, grs_y), (f18_x, f18_y), (np.log10(r21_surf_den.value), np.log10(r21_sigma))

def plot_equilibrium_lines(ax):
    """Plot virial equilibrium and pressure-bounded equilibrium lines."""
    Gamma = 0.73
    G = const.G.value
    k = const.k_B.value
    pc = const.pc.value
    r_range = np.logspace(-2.5, 1.5, num=5000)
    
    # Standard equilibrium
    Pe = 0
    Y = (1/3) * ((np.pi*Gamma*G*r_range*10) + ((4*(Pe*(k*1e6)))/(r_range*10)))
    plt.plot(np.log10(r_range), np.log10((Y*pc)/1e6), color='k', ls='--', alpha=0.2)
    
    # Pressure equilibrium
    pressures = [1e5, 1e6, 1e7, 1e8, 1e9]
    for Pe in pressures:
        Y = (1/3) * ((np.pi*Gamma*G*r_range*10) + ((4*(Pe*(k*1e6)))/(r_range*10)))
        plt.plot(np.log10(r_range), np.log10((Y*pc)/1e6), color='k', alpha=0.2)
    
    plt.text(-1.3, 2.1, r'10$^{8}$ K cm$^{-3}$', color='k', fontsize=10)
    plt.text(-1.7, 1.5, r'10$^{7}$ K cm$^{-3}$', color='k', fontsize=10)
    plt.text(-2.2, 1.0, r'10$^{6}$ K cm$^{-3}$', color='k', fontsize=10)

def create_plot(cmz_data, grs_data, f18_data, r21_data):
    """Create and save the final plot."""
    fig = plt.figure(figsize=(6.5, 5.5))
    
    def scatter_with_outline(x, y, color, label, alpha=0.3, zorder=2):
        plt.scatter(x, y, ec='none', fc='black', s=43, zorder=zorder)
        plt.scatter(x, y, ec='none', fc='white', s=30, zorder=zorder)
        plt.scatter(x, y, c=color, ec='black' if color=='c' else 'none', 
                   s=30, alpha=alpha, label=label, zorder=zorder)
    
    scatter_with_outline(*r21_data, 'green', "(Rosolowsky et al. 2021)", alpha=0.3, zorder=1)
    scatter_with_outline(*cmz_data, 'c', "Milky Way CMZ", alpha=0.6)
    scatter_with_outline(*grs_data, 'red', "Milky Way Clouds", alpha=0.3)
    scatter_with_outline(*f18_data, 'blue', "NGC 300", alpha=0.4)
    
    plot_equilibrium_lines(plt.gca())
    
    plt.ylabel(r'log $\sigma^{2}$/ R (km$^{2}$ s$^{-2}$ pc$^{-1}$)', 
              fontname='Georgia', fontsize=13)
    plt.xlabel(r'log $\Sigma$ (g cm$^{-2}$)', fontname='Georgia', fontsize=13)
    plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=8)
    plt.ylim(-1.0, 3.0)
    plt.xlim(-2.5, 0.5)
    plt.grid(alpha=0.5)
    
    mpl.rc('font', family='Georgia')
    mpl.rc('xtick', labelsize=15)
    mpl.rc('ytick', labelsize=15)
    
    plt.savefig('size_lw_vs_Sigma.pdf')
    plt.clf()

def main():
    # Clouds to drop from the data (not in the CMZ)
    drop_leaves = [0, 1, 2, 3, 4, 5, 6, 7, 10, 55, 56]
    
    kin_data, phys_data = load_and_clean_data(
        'cmz_cloud_kinematic_properties_sorted_final.tex',
        'cmz_cloud_physical_properties_sorted_final.tex',
        drop_leaves
    )
    
    cmz_x, cmz_y = calculate_properties(phys_data, kin_data)
    
    grs_data, f18_data, r21_data = load_comparison_data(
        './data/GRS.txt',
        './data/f18.txt',
        'GMC_catalogue.ecsv'
    )
    
    create_plot((cmz_x, cmz_y), grs_data, f18_data, r21_data)

if __name__ == "__main__":
    main()

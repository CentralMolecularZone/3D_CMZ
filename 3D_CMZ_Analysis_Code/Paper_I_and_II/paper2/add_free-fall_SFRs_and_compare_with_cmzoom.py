import math
import pandas as pd
import numpy as np
from astropy.table import Table
from scipy.integrate import quad
from astropy import units as u
from astropy import constants as const


def IMF_powerlaw_mass(M, alpha, k):
	M_tot = (k * M * M**(-alpha))
	return M_tot


physical_table = Table.read('dendro_physical.tex')
mass = (physical_table['Mass']*u.M_sun).to(u.kg)
radius = (physical_table['Radius']*u.pc).to(u.m)
rho = mass/((4/3)*np.pi*(radius**3))
tff = np.sqrt(3*np.pi/(32*const.G*rho))
luminosity_table = Table.read('dendro_sfr.tex')
ltot = luminosity_table['L$_{\\textrm{total IR}}$']

# --- # IMF definitions
mass_star_tot_temp = np.empty(3) * np.nan
M_upper_lim = np.inf
N=1.0
alpha1 = 2.3
alpha2 = 1.3
alpha3 = 0.3

# Add new SFRs based on Barnes et al. 2017, using free-fall time as the time scale
SFR =[]
for i in range(len(tff)):
	sum_lum_log = math.log10(ltot[i])
	sum_mass_star = (sum_lum_log**3 * 4.77751404) + (sum_lum_log**2 * -52.637041) + (sum_lum_log*198.602281) - 245.57466
	if sum_mass_star > 0:
		k1 = (N*1.3) * (sum_mass_star**(-1.3) - (M_upper_lim)**(-1.3))**-1
		mass_star_tot_temp[2] = quad(IMF_powerlaw_mass, 0.5, M_upper_lim, args=(alpha1,k1))[0]
		k2 = k1 * (0.5**(-alpha1))/(0.5**(-alpha2))
		mass_star_tot_temp[1] = quad(IMF_powerlaw_mass, 0.08, 0.5, args=(alpha2,k2))[0]
		k3 = k2 * (0.08**(-alpha2))/0.08**(-alpha3)
		mass_star_tot_temp[0] = quad(IMF_powerlaw_mass, 0.001, 0.08, args=(alpha3,k3))[0]
		mass_star_tot = np.nansum(mass_star_tot_temp)
		SFR.append(mass_star_tot/(tff[i].to(u.yr).value))
	else:
		SFR.append(np.nan)

luminosity_table['SFR_Barnes17'] = SFR
luminosity_table['SFR_Barnes17'].format = '%.1E'
luminosity_table['SFR$_{\\textrm{total IR}}$'].format = '%.1E'
luminosity_table.write('dendro_sfr.tex', format='ascii.latex', overwrite=True)


# Create a table comparing the SFRs from this paper to the SFRs from the CMZoom survey
cmzoom_tab = pd.read_csv('Herschel_CMZ_CMZoom_SFR_comparison.csv')
sf_tab = luminosity_table.to_pandas()
sf_tab.rename(columns={'SFR$_{\\textrm{total IR}}$': 'SFR_TotalIR'}, inplace=True)
sf_tab['SFR_TotalIR'] = sf_tab['SFR_TotalIR']*1e3
sf_tab['SFR_Barnes17'] = sf_tab['SFR_Barnes17']*1e3
sf_tab_comparison = sf_tab[['Leaf ID', 'SFR_TotalIR', 'SFR_Barnes17']]

cols_to_append = ['Robust CMZoom SFR (Herschel Temp)', 'Robust CMZoom SFR (50K)', 
                  'Robust and ambig CMZoom SFR (Herschel Temp)', 'Robust and ambig CMZoom SFR (50K)']

result = sf_tab_comparison.merge(cmzoom_tab[['Herschel Mask'] + cols_to_append], 
                                 left_on='Leaf ID', 
                                 right_on='Herschel Mask', 
                                 how='left').drop(columns='Herschel Mask')

result['SFR_TotalIR'] = result['SFR_TotalIR'].round(2)
result['SFR_Barnes17'] = result['SFR_Barnes17'].round(2)
result.to_csv('sf_tab_comparison.csv', index=False)
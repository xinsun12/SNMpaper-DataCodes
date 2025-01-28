#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:24:07 2024

@author: xinsun
"""

#%% imports
import sys
import os
import numpy as np
import xarray as xr
import pandas as pd

# plotting packages
import seaborn as sb
sb.set(style='ticks')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.gridspec import GridSpec
import cmocean.cm as cmo
from cmocean.tools import lighten

# numerical packages
from numba import jit

#%% read data from Virtual chemostat
os.chdir("/yourpath/output0D") 
#os.chdir("/yourpath/output0D_withoutNOB") 

fname = 'OMpulse_P0.16'

out_O2 = np.loadtxt(fname+'_O2.txt', delimiter='\t')
out_NO2 = np.loadtxt(fname+'_NO2.txt', delimiter='\t')
out_Sd = np.loadtxt(fname+'_OM.txt', delimiter='\t')
out_bHet = np.loadtxt(fname+'_bHet.txt', delimiter='\t')
out_b1Den = np.loadtxt(fname+'_b1Den.txt', delimiter='\t')
out_b2Den = np.loadtxt(fname+'_b2Den.txt', delimiter='\t')
out_b4Den = np.loadtxt(fname+'_b4Den.txt', delimiter='\t')
out_bAOO = np.loadtxt(fname+'_bAOO.txt', delimiter='\t')
out_bNOO = np.loadtxt(fname+'_bNOO.txt', delimiter='\t')
out_bAOX = np.loadtxt(fname+'_bAOX.txt', delimiter='\t')


meanNO2 = np.loadtxt(fname+'_meanNO2.txt', delimiter='\t')

days = 4e4
dt = 0.001
timesteps = days/dt
out_at_day = dt
nn_output = days/out_at_day
nn_outputforaverage = int(2000/out_at_day) 

#%% Figure 1
fstic = 13
fslab = 15
fig2 = plt.figure(figsize=(4,5))
gs = GridSpec(2,1)

lineplotx = np.arange(nn_output-100/out_at_day, nn_output+1, 1)*out_at_day-39900
                    
ax1 = plt.subplot(gs[0,0])
plt.title('Nutrient (µM)', fontsize = fslab)
plt.plot(lineplotx, out_Sd, color='firebrick', linestyle='-', linewidth = 2, label='OM')
plt.plot(lineplotx, out_NO2, color='royalblue', linestyle='-', linewidth = 2, label='NO$_2$')
plt.axhline(y = meanNO2, color='royalblue', linestyle='--', linewidth = 3, label='mean NO$_2$')  
plt.plot(lineplotx, out_O2, color='grey', linestyle='-', linewidth = 2, label='O$_2$')             
ax1.set_ylim([0.0, 8])
ax1.set_yticks(np.arange(0.0,9,2))   

           
ax2 = plt.subplot(gs[1,0])
plt.title('Biomass (µM-N)', fontsize = fslab)
plt.plot(lineplotx, out_bHet, color='grey', linestyle='--', linewidth = 2,label='aer')
plt.plot(lineplotx, out_b1Den, color='firebrick', linestyle='-', linewidth = 2,label='NO$_3$→NO$_2$')
plt.plot(lineplotx, out_b4Den, color='goldenrod', linestyle='-', linewidth = 2,label='NO$_2$→N$_2$O')
plt.plot(lineplotx, out_b2Den, color='goldenrod', linestyle='--', linewidth = 2,label='NO$_2$→N$_2$')
plt.plot(lineplotx, out_bAOO, color='grey', linestyle='-', linewidth = 2,label='AOO')
plt.plot(lineplotx, out_bNOO, color='royalblue', linestyle='-', linewidth = 2,label='NOO')
plt.plot(lineplotx, out_bAOX, color='forestgreen', linestyle='-', linewidth = 2,label='AOX')

ax2.set_ylim([0.0, 0.3])
ax2.set_yticks(np.arange(0.0,0.31,0.1))
                
## delete axis title of some subplots
ax1.tick_params(labelsize=fstic, labelbottom=False)
ax2.tick_params(labelsize=fstic)

ax2.set_xlabel('days', fontsize=fslab)

# select part of x axis 
xlowerlimit = 0
xupperlimit = 100
xtickdiv = 20
for ax in [ax1, ax2]:
    ax.set_xlim([xlowerlimit, xupperlimit])
    ax.set_xticks(np.arange(xlowerlimit, xupperlimit+1, xtickdiv))  
                
plt.tight_layout()

#%% Save plots             
os.chdir("/yourpath")
plt.savefig('Fig1.png', dpi=300)
#plt.savefig('Fig1_noNOB.png', dpi=300)

#%% read data from ROMS results
# Load the data from Excel
file_path = '/yourpath/nitrox_20years.xlsx'  # Update with the path to your Excel file
# # Notes
# DENITRIF1 = NO3 --> NO2 (via NAR)
# DENITRIF2 = NO2 --> N2O (via NIR)
# DENITRIF3 = N2O --> N2 (via NOS)
# DENITRIF4 = NO3 --> N2O (via NAI)
# DENITRIF5 = NO2 --> N2 (via NIO)
# DENITRIF6 = NO3 --> N2 (via NAO)

o2_mean = pd.read_excel(file_path, sheet_name='O2_mean', index_col=0)
o2_last_5_years = o2_mean.iloc[:, -60:]
o2_onemean = o2_last_5_years.mean(axis=1)

no2_mean = pd.read_excel(file_path, sheet_name='NO2_mean', index_col=0)
no2_last_5_years = no2_mean.iloc[:, -60:]
no2_onemean = no2_last_5_years.mean(axis=1)


no2n2oR_mean = pd.read_excel(file_path, sheet_name='DENITRIF2_mean', index_col=0)
no2n2oR_last_5_years = no2n2oR_mean.iloc[:, -60:]
no2n2oR_onemean = no2n2oR_last_5_years.mean(axis=1)

no2n2R_mean = pd.read_excel(file_path, sheet_name='DENITRIF5_mean', index_col=0)
no2n2R_last_5_years = no2n2R_mean.iloc[:, -60:]
no2n2R_onemean = no2n2R_last_5_years.mean(axis=1)

no3no2R_mean = pd.read_excel(file_path, sheet_name='DENITRIF1_mean', index_col=0)
no3no2R_last_5_years = no3no2R_mean.iloc[:, -60:]
no3no2R_onemean = no3no2R_last_5_years.mean(axis=1)

nobR_mean = pd.read_excel(file_path, sheet_name='NITROX_mean', index_col=0)
nobR_last_5_years = nobR_mean.iloc[:, -60:]
nobR_onemean = nobR_last_5_years.mean(axis=1)

aoaR_mean = pd.read_excel(file_path, sheet_name='AMMOX_mean', index_col=0)
aoaR_last_5_years = aoaR_mean.iloc[:, -60:]
aoaR_onemean = aoaR_last_5_years.mean(axis=1)

aoxR_mean = pd.read_excel(file_path, sheet_name='ANAMMOX_mean', index_col=0)
aoxR_last_5_years = aoxR_mean.iloc[:, -60:]
aoxR_onemean = aoxR_last_5_years.mean(axis=1)

#%% read data from observations
nitrite_data = pd.read_csv('/yourpath/ObsConcentrations.csv')
rates_data = pd.read_csv('/yourpath/ObsRates.csv')
#%% Function for plot with shaded areas
def plot_with_shading(ax, data, x_mean, x_std, label, color, linestyle='-', flip_sign=False):
    """
    Plots a line with shaded areas representing standard deviations.

    Parameters:
    - data: Index values (e.g., depth or time) for the y-axis
    - x_mean: Mean values for the x-axis
    - x_std: Standard deviation values for the shaded area
    """
    x_values = -x_mean if flip_sign else x_mean
    ax.plot(x_values, data, linestyle, color=color, linewidth=2, label=label)
    ax.fill_betweenx(data, 
                     x_values - x_std, 
                     x_values + x_std, 
                     color=color, alpha=0.3, linewidth=0)


#%% Plot depth profiles of ROMs and observations with shaded areas
# Plot
plt.figure(figsize=(6, 6))
gs = GridSpec(2, 2)
shade_alpha = 0.3
point_alpha = 0.8
pointsize = 20
fstic = 13
fslab = 15 
# Plot 1
ax1 = plt.subplot(gs[0, 0])
plt.plot(no2_onemean*100, no2_last_5_years.index, '-', color='royalblue', linewidth=2, label='NO$_2$$^–$x100')
plt.fill_betweenx(no2_last_5_years.index, 
                  (no2_onemean - 2*no2_last_5_years.std(axis=1)) * 100,
                  (no2_onemean + 2*no2_last_5_years.std(axis=1)) * 100,
                  color='royalblue', linewidth=0, alpha=0.3)
plt.plot(o2_onemean, o2_last_5_years.index, '-', color='grey', linewidth=2, label='O$_2$')
plt.fill_betweenx(o2_last_5_years.index, 
                  o2_onemean - 2*o2_last_5_years.std(axis=1),
                  o2_onemean + 2*o2_last_5_years.std(axis=1),
                  color='grey',linewidth=0, alpha=0.3)


plt.ylim([0.0, 600])
plt.yticks(np.arange(0.0, 601, 100)) 
plt.xlim([0.0, 380])
plt.xticks(np.arange(0.0, 301, 100)) 
plt.gca().invert_yaxis()
#plt.legend()
plt.ylabel('Depth (m)', fontsize=fslab)

# Plot 2
ax2 = plt.subplot(gs[0, 1])

# Vertical line at x=0
ax2.axvline(x=0, color='k', linestyle='-', linewidth=1)
plot_with_shading(ax2, nobR_last_5_years.index, -nobR_onemean, 2*nobR_last_5_years.std(axis=1), 'NOB', 'royalblue') #nobR_onestd, 
plot_with_shading(ax2, no2n2R_last_5_years.index, -no2n2R_onemean, 2*no2n2R_last_5_years.std(axis=1), 'NO2-->N2', 'goldenrod', linestyle='--')
plot_with_shading(ax2, no2n2oR_last_5_years.index, -no2n2oR_onemean, 2*no2n2oR_last_5_years.std(axis=1), 'NO2-->N2O', 'goldenrod')
plot_with_shading(ax2, no3no2R_last_5_years.index, no3no2R_onemean, 2*no3no2R_last_5_years.std(axis=1), 'NO3-->NO2', 'firebrick')
plot_with_shading(ax2, aoaR_last_5_years.index, aoaR_onemean, 2*aoaR_last_5_years.std(axis=1), 'AOA', 'grey')
plot_with_shading(ax2, aoxR_last_5_years.index, -aoxR_onemean, 2*aoxR_last_5_years.std(axis=1), 'AOX', 'forestgreen')

plt.xlim([-0.16, 0.16])
plt.ylim([0.0, 600])
plt.yticks(np.arange(0.0, 601, 100)) 
plt.gca().invert_yaxis() 

# Plot 3
ax3 = plt.subplot(gs[1, 0])
for cast in nitrite_data['Cast'].unique():
    cast_data = nitrite_data[nitrite_data['Cast'] == cast]
    ax3.plot(cast_data['Nitrite'] * 100, cast_data['Depth'], '-',
             color="royalblue", alpha=1, linewidth=2)
    ax3.scatter(cast_data['Nitrite'] * 100, cast_data['Depth'], 
                color="royalblue", s=pointsize, alpha=point_alpha)
for cast in nitrite_data['Cast'].unique():
    cast_data = nitrite_data[nitrite_data['Cast'] == cast]
    ax3.plot(cast_data['Sbeox0Mm/L'], cast_data['Depth'], '-',
             color="grey", alpha=1, linewidth=2)
    ax3.scatter(cast_data['Sbeox0Mm/L'], cast_data['Depth'], 
                color="grey", s=pointsize, alpha=point_alpha)
plt.ylim([0.0, 600])
plt.yticks(np.arange(0.0, 601, 100)) 
plt.xlim([0.0, 380])
plt.xticks(np.arange(0.0, 301, 100)) 
plt.gca().invert_yaxis()
plt.ylabel('Depth (m)', fontsize=fslab)
plt.xlabel('Concentration (µM)', fontsize=fslab)

# Plot 4
ax4 = plt.subplot(gs[1, 1])
ax4.axvline(x=0, color='k', linestyle='-', linewidth=1)
ax4.fill_betweenx(rates_data['Depth_m'], 
                  (rates_data['NO3_reduction_nM_N/d'] - rates_data['Error.1']) / 1000,
                  (rates_data['NO3_reduction_nM_N/d'] + rates_data['Error.1']) / 1000,
                  color="firebrick", linewidth=0, alpha=shade_alpha)
ax4.plot(rates_data['NO3_reduction_nM_N/d'] / 1000, rates_data['Depth_m'], '-',
         color="firebrick", linewidth=2)
ax4.scatter(rates_data['NO3_reduction_nM_N/d'] / 1000, rates_data['Depth_m'], 
            color="firebrick", s=pointsize, alpha=point_alpha)
# no2-->n2 rates
ax4.fill_betweenx(rates_data['Depth_m'], 
                  -(rates_data['DN_Rate_nM_N2/day'] - rates_data['Error.4']) * 2 / 1000,
                  -(rates_data['DN_Rate_nM_N2/day'] + rates_data['Error.4']) * 2 / 1000,
                  color="goldenrod", linewidth=0, alpha=shade_alpha)
ax4.plot(-(rates_data['DN_Rate_nM_N2/day'] * 2 / 1000), rates_data['Depth_m'], '--',
         color="goldenrod", linewidth=2)
ax4.scatter(-(rates_data['DN_Rate_nM_N2/day'] * 2 / 1000), rates_data['Depth_m'], 
            color="goldenrod", s=pointsize, alpha=point_alpha)
#Anammox
ax4.fill_betweenx(rates_data['Depth_m'], 
                  -(rates_data['AMX_Rate_nM_N2/day'] - rates_data['Error.3']) * 2 / 1000,
                  -(rates_data['AMX_Rate_nM_N2/day'] + rates_data['Error.3']) * 2 / 1000,
                  color="forestgreen", linewidth=0, alpha=shade_alpha)
ax4.plot(-(rates_data['AMX_Rate_nM_N2/day'] * 2 / 1000), rates_data['Depth_m'], '-',
         color="forestgreen", linewidth=2)
ax4.scatter(-(rates_data['AMX_Rate_nM_N2/day'] * 2 / 1000), rates_data['Depth_m'], 
            color="forestgreen", s=pointsize, alpha=point_alpha)
# Nitrite oxi
ax4.fill_betweenx(rates_data['Depth_m'], 
                  -(rates_data['Xin_NO2_oxi_rate_nM/d'] - rates_data['se']) / 1000,
                  -(rates_data['Xin_NO2_oxi_rate_nM/d'] + rates_data['se']) / 1000,
                  color="royalblue", linewidth=0, alpha=shade_alpha)
ax4.plot(-(rates_data['Xin_NO2_oxi_rate_nM/d'] / 1000), rates_data['Depth_m'], '-',
         color="royalblue", linewidth=2)
ax4.scatter(-(rates_data['Xin_NO2_oxi_rate_nM/d'] / 1000), rates_data['Depth_m'], 
            color="royalblue", s=pointsize, alpha=point_alpha)
# Ammonia oxi
ax4.fill_betweenx(rates_data['Depth_m'], 
                  (rates_data['NH4_oxi_nM_N/d'] - rates_data['Error']) / 1000,
                  (rates_data['NH4_oxi_nM_N/d'] + rates_data['Error']) / 1000,
                  color="grey", linewidth=0, alpha=shade_alpha)
ax4.plot(rates_data['NH4_oxi_nM_N/d'] / 1000, rates_data['Depth_m'], '-',
         color="grey", linewidth=2)
ax4.scatter(rates_data['NH4_oxi_nM_N/d'] / 1000, rates_data['Depth_m'], 
            color="grey", s=pointsize, alpha=point_alpha)

plt.xlim([-0.16, 0.16])
plt.ylim([0.0, 600])
plt.yticks(np.arange(0.0, 601, 100)) 
plt.gca().invert_yaxis() 
plt.xlabel('Rate (µM d$^-$$^1$)', fontsize=fslab)


ax1.tick_params(labelsize=fstic, labelbottom=False)
ax2.tick_params(labelsize=fstic, labelbottom=False, labelleft=False)
ax4.tick_params(labelsize=fstic)
ax4.tick_params(labelsize=fstic, labelleft=False)
plt.tight_layout()
plt.show()

#%% Save plots             
os.chdir("/yourpath")
plt.savefig('Fig2_depthprofiles.png', dpi=300)

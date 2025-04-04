# -*- coding: utf-8 -*-
"""
Created in March 2023, last update: April 2025

Purpose
-------
    A 0D chemostat model with microbes in marine OMZs with a focus on nitrite accumulation,
    Modular denitrification included, yields of denitrifiers depend on Gibbs free energy.
    Organic matter pulses included
    
@authors: 
    Xin Sun, Emily Zakem, Pearse Buchanan
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

#%% Set initial conditions and incoming concentrations to chemostat experiment

### Organic matter (S)
## Constant supply (assuming 1 m3 box and flux into top of box)
Sd0_exp = 0.05 #µM-N m-3 day-1
## Pulse intensity 
xpulse_Sd = np.array([2]) #µM-N

### Oxygen supply
O20_exp =  np.array([2])

### model parameters for running experiments
### dil = 0 for checking N balance 
dil = 0.04  # dilution rate (1/day)
if dil == 0:
    days = 10  # number of days to run chemostat
    dt = 0.001  # timesteps per day (days)
    timesteps = days/dt     # number of timesteps
    out_at_day = 0.1       # output results this often (days)
    nn_output = days/out_at_day     # number of entries for output
    print("dilution = 0, check N balance")
else:
    days = 5e4  # number of days to run chemostat
    dt = 0.001  # timesteps length (days)
    timesteps = days/dt     # number of timesteps
    out_at_day = dt         # output results this often
    nn_output = days/out_at_day     # number of entries for output
    print("dilution > 0, run experiments")
    
nn_outputforaverage = int(2000/out_at_day) # finish value is the average of the last XX (number) of outputs
     
#%% Define variables  
outputd1 = xpulse_Sd 
outputd2 = O20_exp

#%% initialize arrays for output

# Nutrients
fin_O2 = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_Sd = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_NO3 = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_NO2 = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_NH4 = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_N2 = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_N2O = np.ones((len(outputd1), len(outputd2))) * np.nan 
# Biomasses
fin_bHet = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_b1Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_b2Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_b3Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_b4Den = np.ones((len(outputd1), len(outputd2))) * np.nan 
fin_b5Den = np.ones((len(outputd1), len(outputd2))) * np.nan 
fin_b6Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_b7Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_bAOO = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_bNOO = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_bAOX = np.ones((len(outputd1), len(outputd2))) * np.nan
# Growth rates
fin_uHet = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_u1Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_u2Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_u3Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_u4Den = np.ones((len(outputd1), len(outputd2))) * np.nan 
fin_u5Den = np.ones((len(outputd1), len(outputd2))) * np.nan 
fin_u6Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_u7Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_uAOO = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_uNOO = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_uAOX = np.ones((len(outputd1), len(outputd2))) * np.nan
# Rates 
fin_rHet = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_rHetAer = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_rO2C = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_r1Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_r2Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_r3Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_r4Den = np.ones((len(outputd1), len(outputd2))) * np.nan 
fin_r5Den = np.ones((len(outputd1), len(outputd2))) * np.nan 
fin_r6Den = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_rAOO = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_rNOO = np.ones((len(outputd1), len(outputd2))) * np.nan
fin_rAOX = np.ones((len(outputd1), len(outputd2))) * np.nan

#%% set traits of the different biomasses
os.chdir("yourpath/ChemostatModel_ModularDenitrification_clean")
### set traits
from traits import * 

#%% calculate R*-stars
from O2_star_Xin import O2_star
from N2O_star_Xin import N2O_star 
from R_star_Xin import R_star

# O2 (nM-O2) 
O2_star_aer = R_star(dil, K_o2_aer, mumax_Het / y_oO2, y_oO2) * 1e3 
O2_star_aoo = R_star(dil, K_o2_aoo, mumax_AOO / y_oAOO, y_oAOO) * 1e3 
O2_star_noo = R_star(dil, K_o2_noo, mumax_NOO / y_oNOO, y_oNOO) * 1e3 
# N2O (nM-N)
N2O_star_den5 = R_star(dil, K_n2o_Den, VmaxN_5Den, y_n5N2O) * 1e3
# OM
OM_star_aer = R_star(dil, K_s, VmaxS, y_oHet)
OM_star_den1 = R_star(dil, K_s, VmaxS, y_n1Den)
OM_star_den2 = R_star(dil, K_s, VmaxS, y_n2Den)
OM_star_den3 = R_star(dil, K_s, VmaxS, y_n3Den)
OM_star_den4 = R_star(dil, K_s, VmaxS, y_n4Den)
OM_star_den5 = R_star(dil, K_s, VmaxS, y_n5Den)
OM_star_den6 = R_star(dil, K_s, VmaxS, y_n6Den)
# Ammonia
Amm_star_aoo = R_star(dil, K_n_AOO, VmaxN_AOO, y_nAOO)
Amm_star_aox = R_star(dil, K_nh4_AOX, VmaxNH4_AOX, y_nh4AOX)
# Nitrite
nitrite_star_den2 = R_star(dil, K_n_Den, VmaxN_2Den, y_n2NO2)
nitrite_star_den4 = R_star(dil, K_n_Den, VmaxN_4Den, y_n4NO2)
nitrite_star_noo = R_star(dil, K_n_NOO, VmaxN_NOO, y_nNOO)
nitrite_star_aox = R_star(dil, K_no2_AOX, VmaxNO2_AOX, y_no2AOX)
# Nitrate
nitrate_star_den1 = R_star(dil, K_n_Den, VmaxN_1Den, y_n1NO3)
nitrate_star_den3 = R_star(dil, K_n_Den, VmaxN_3Den, y_n3NO3)
nitrate_star_den6 = R_star(dil, K_n_Den, VmaxN_6Den, y_n6NO3)


#%% begin loop of experiments

from model import OMZredox

for k in np.arange(len(outputd1)):
    for m in np.arange(len(O20_exp)):
        print(k,m)
        
        # 1) Chemostat influxes (µM-N or µM O2)
        in_Sd = Sd0_exp
        in_O2 = O20_exp[m]
        in_NO3 = 30.0
        in_NO2 = 0.0
        in_NH4 = 0.0
        in_N2 = 0.0
        in_N2O = 0.0
        # initial conditions
        initialOM = in_Sd
        initialNO2 = 0
        
        # 2) Initial biomasses (set to 0.0 to exclude a microbial group, 0.1 as default) 
        in_bHet = 0.1
        in_b1Den = 0.1 # NO3-->NO2, cross-feed
        in_b4Den = 0.1 # NO2-->N2O, cross-feed
        in_b5Den = 0.1 # N2O-->N2, cross-feed
        in_b2Den = 0.1 # NO2-->N2
        in_b3Den = 0.1 # complete denitrifier
        in_b6Den = 0.1 # NO3-->N2O
        in_b7Den = 0 # bookend: NO3-->NO2, N2O-->N2    
        in_bAOO = 0.1
        in_bNOO = 0.1
        in_bAOX = 0.1
        
        # pulse conditions        
        pulse_int = 50 #pulse interval
        pulse_Sd = xpulse_Sd[k]#pulse intensity
        pulse_O2 = 0.0
       
       
        # 3) Call main model
        results = OMZredox(timesteps, nn_output, dt, dil, out_at_day, \
                           pulse_Sd, pulse_O2, pulse_int, \
                           K_o2_aer, K_o2_aoo, K_o2_noo, \
                           K_n2o_Den, \
                           mumax_Het, mumax_AOO, mumax_NOO, mumax_AOX, \
                           VmaxS, K_s, \
                           VmaxN_1Den, VmaxN_2Den, VmaxN_3Den, VmaxN_4Den, VmaxN_5Den, VmaxN_6Den, K_n_Den, \
                           VmaxN_AOO, K_n_AOO, VmaxN_NOO, K_n_NOO, \
                           VmaxNH4_AOX, K_nh4_AOX, VmaxNO2_AOX, K_no2_AOX, \
                           y_oHet, y_oO2, \
                           y_n1Den, y_n1NO3, y_n2Den, y_n2NO2, y_n3Den, y_n3NO3, y_n4Den, y_n4NO2, y_n5Den, y_n5N2O, y_n6Den, y_n6NO3, y_n7Den_NO3, y_n7NO3, e_n7Den_NO3, y_n7Den_N2O, y_n7N2O, e_n7Den_N2O,\
                           y_nAOO, y_oAOO, y_nNOO, y_oNOO, y_nh4AOX, y_no2AOX, \
                           e_n2Den, e_n3Den, e_no3AOX, e_n2AOX, e_n4Den, e_n5Den, e_n6Den, e_n1Den, \
                           initialOM, initialNO2, in_Sd, in_O2, in_NO3, in_NO2, in_NH4, in_N2, in_N2O, \
                           in_bHet, in_b1Den, in_b2Den, in_b3Den, in_bAOO, in_bNOO, in_bAOX, in_b4Den, in_b5Den, in_b6Den, in_b7Den)
        
        out_Sd = results[0]
        out_O2 = results[1]
        out_NO3 = results[2]
        out_NO2 = results[3]
        out_NH4 = results[4]
        out_N2O = results[5] 
        out_N2 = results[6]
        out_bHet = results[7]
        out_b1Den = results[8]
        out_b2Den = results[9]
        out_b3Den = results[10]
        out_b4Den = results[11]
        out_b5Den = results[12]
        out_b6Den = results[13]
        out_bAOO = results[14]
        out_bNOO = results[15]
        out_bAOX = results[16]
        out_uHet = results[17]
        out_u1Den = results[18]
        out_u2Den = results[19]
        out_u3Den = results[20]
        out_u4Den = results[21]
        out_u5Den = results[22]
        out_u6Den = results[23]
        out_uAOO = results[24]
        out_uNOO = results[25]
        out_uAOX = results[26]
        out_rHet = results[27]
        out_rHetAer = results[28]
        out_rO2C = results[29]
        out_r1Den = results[30]
        out_r2Den = results[31]
        out_r3Den = results[32]
        out_r4Den = results[33]
        out_r5Den = results[34]
        out_r6Den = results[35]
        out_rAOO = results[36]
        out_rNOO = results[37]
        out_rAOX = results[38]     
        out_b7Den = results[39]
        out_u7Den = results[40]
        
        ## 4) plot line plots
        if len(outputd1)*len(outputd2) < 3:
            if dil == 0:
                print("test run for dilution rate = 0")
            else:
                fin_NO2[k,m] = np.nanmean(out_NO2[-nn_outputforaverage::])
                os.chdir("/yourpath/output0D") 

                fname = 'OMpulse_P0.16'
                np.savetxt(fname+'_meanNO2.txt', fin_NO2[0], delimiter='\t')
                np.savetxt(fname+'_O2.txt', out_O2[-100001::], delimiter='\t')
                np.savetxt(fname+'_NO2.txt', out_NO2[-100001::], delimiter='\t')
                np.savetxt(fname+'_OM.txt', out_Sd[-100001::], delimiter='\t')

                np.savetxt(fname+'_bHet.txt', out_bHet[-100001::], delimiter='\t')
                np.savetxt(fname+'_b1Den.txt', out_b1Den[-100001::], delimiter='\t')
                np.savetxt(fname+'_b2Den.txt', out_b2Den[-100001::], delimiter='\t')
                np.savetxt(fname+'_b4Den.txt', out_b4Den[-100001::], delimiter='\t')
                np.savetxt(fname+'_bAOO.txt', out_bAOO[-100001::], delimiter='\t')
                np.savetxt(fname+'_bNOO.txt', out_bNOO[-100001::], delimiter='\t')
                np.savetxt(fname+'_bAOX.txt', out_bAOX[-100001::], delimiter='\t')
        else:
            print("experiments >= 3, turn off lineplots")
    
        # 5) Record solutions in initialised arrays
        fin_O2[k,m] = np.nanmean(out_O2[-nn_outputforaverage::])
        fin_Sd[k,m] = np.nanmean(out_Sd[-nn_outputforaverage::])
        fin_NO3[k,m] = np.nanmean(out_NO3[-nn_outputforaverage::])
        fin_NO2[k,m] = np.nanmean(out_NO2[-nn_outputforaverage::])
        fin_NH4[k,m] = np.nanmean(out_NH4[-nn_outputforaverage::])
        fin_N2[k,m] = np.nanmean(out_N2[-nn_outputforaverage::])
        fin_N2O[k,m] = np.nanmean(out_N2O[-nn_outputforaverage::]) 
        fin_bHet[k,m] = np.nanmean(out_bHet[-nn_outputforaverage::])
        fin_b1Den[k,m] = np.nanmean(out_b1Den[-nn_outputforaverage::])
        fin_b2Den[k,m] = np.nanmean(out_b2Den[-nn_outputforaverage::])
        fin_b3Den[k,m] = np.nanmean(out_b3Den[-nn_outputforaverage::])
        fin_b4Den[k,m] = np.nanmean(out_b4Den[-nn_outputforaverage::]) 
        fin_b5Den[k,m] = np.nanmean(out_b5Den[-nn_outputforaverage::]) 
        fin_b6Den[k,m] = np.nanmean(out_b6Den[-nn_outputforaverage::])
        fin_b7Den[k,m] = np.nanmean(out_b7Den[-nn_outputforaverage::]) 
        fin_bAOO[k,m] = np.nanmean(out_bAOO[-nn_outputforaverage::])
        fin_bNOO[k,m] = np.nanmean(out_bNOO[-nn_outputforaverage::])
        fin_bAOX[k,m] = np.nanmean(out_bAOX[-nn_outputforaverage::])
        fin_uHet[k,m] = np.nanmean(out_uHet[-nn_outputforaverage::])
        fin_u1Den[k,m] = np.nanmean(out_u1Den[-nn_outputforaverage::])
        fin_u2Den[k,m] = np.nanmean(out_u2Den[-nn_outputforaverage::])
        fin_u3Den[k,m] = np.nanmean(out_u3Den[-nn_outputforaverage::])
        fin_u4Den[k,m] = np.nanmean(out_u4Den[-nn_outputforaverage::]) 
        fin_u5Den[k,m] = np.nanmean(out_u5Den[-nn_outputforaverage::]) 
        fin_u6Den[k,m] = np.nanmean(out_u6Den[-nn_outputforaverage::]) 
        fin_u7Den[k,m] = np.nanmean(out_u7Den[-nn_outputforaverage::])
        fin_uAOO[k,m] = np.nanmean(out_uAOO[-nn_outputforaverage::])
        fin_uNOO[k,m] = np.nanmean(out_uNOO[-nn_outputforaverage::])
        fin_uAOX[k,m] = np.nanmean(out_uAOX[-nn_outputforaverage::])
        fin_rHet[k,m] = np.nanmean(out_rHet[-nn_outputforaverage::])
        fin_rHetAer[k,m] = np.nanmean(out_rHetAer[-nn_outputforaverage::])
        fin_rO2C[k,m] = np.nanmean(out_rO2C[-nn_outputforaverage::])
        fin_r1Den[k,m] = np.nanmean(out_r1Den[-nn_outputforaverage::])
        fin_r2Den[k,m] = np.nanmean(out_r2Den[-nn_outputforaverage::])
        fin_r3Den[k,m] = np.nanmean(out_r3Den[-nn_outputforaverage::])
        fin_r4Den[k,m] = np.nanmean(out_r4Den[-nn_outputforaverage::]) 
        fin_r5Den[k,m] = np.nanmean(out_r5Den[-nn_outputforaverage::]) 
        fin_r6Den[k,m] = np.nanmean(out_r6Den[-nn_outputforaverage::]) 
        fin_rAOO[k,m] = np.nanmean(out_rAOO[-nn_outputforaverage::])
        fin_rNOO[k,m] = np.nanmean(out_rNOO[-nn_outputforaverage::])
        fin_rAOX[k,m] = np.nanmean(out_rAOX[-nn_outputforaverage::])
    
            
# delete results only save fin (average)
del results
del out_Sd, out_O2, out_NO3, out_NO2, out_NH4, out_N2, out_N2O
del out_bHet, out_b1Den, out_b2Den, out_b3Den, out_b4Den, out_b5Den, out_b6Den, out_b7Den, out_bAOO, out_bNOO, out_bAOX
del out_uHet, out_u1Den, out_u2Den, out_u3Den, out_u4Den, out_u5Den, out_u6Den, out_u7Den, out_uAOO, out_uNOO, out_uAOX
del out_rHet, out_rHetAer, out_rO2C, out_r1Den, out_r2Den, out_r3Den, out_r4Den, out_r5Den, out_r6Den, out_rAOO, out_rNOO, out_rAOX




## This repository contains Python, Matlab, and Fortran scripts and data used to generate the figures presented in the paper "Microbial Interactions Matter for Biogeochemistry: Nitrite Consumers Induce Nitrite Accumulation"

### 1. Scripts
#### A. Scripts for model runs (outputs of these model runs are available, see 2. Data files below)
- Folder '0D_files':
  All Python scripts containing the code used to run the 0-D model.
- Folder 'ROMS_files':
  All Fortran scripts containing the code used to run ROMS.
#### B. Scripts for making figures
- SNM_Fig1allFig2profiles.py:
  This Python script contains the code used to generate Figure 1a-d (presenting results from the virtual chemostat model), and Figure 4a-d (presenting depth profiles from ROMS and observations).
- Sun_et_al_ROMS_figures.m:
  This Matlab script contains the MATLAB code used to generate Figure 2e&f
  with additional functions in 'scripts/'
  

### 2. Data files
- output0D (results from the virtual chemostat model)
- output0D_withoutNOB (results from the virtual chemostat model excluding NOB)
- nitrox_20years (results from ROMS for Figure 4a&c)
- ROMS output for Figure 2b&d can be found in 'data/'
- ObsConcentrations (measured concentrations from published studies: Tracey et al., 2023; Sun et al., 2023)
- ObsRates (measured concentrations from published studies: Tracey et al., 2023; Sun et al., 2023)

# DALEC-GRASS
This fortran code is a development of the [Data Assimilation Linked Ecosystem Carbon (DALEC)](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/1051-0761%281997%29007%5B0882%3APGPPIT%5D2.0.CO%3B2) model that simulates carbon (C) dynamics in managed grassland ecosystems. DALEC-Grass requires 16 input variables and uses 37 parameters to estimate various C pools and fluxes. Model inputs, parameters, pools and fluxes are described in the first lines of the [code](https://github.com/vmyrgiotis/DALEC_Grass/blob/master/DALEC_GRASS.f90). The [Aggregated Canopy Model (ACM)](https://doi.org/10.1890/1051-0761(1997)007[0882:PGPPIT]2.0.CO;2) lies at the core of DALEC and is used to estimate daily Gross Primary Production (GPP). At each timestep, DALEC allocates the estimated Net Primary Production (NPP) to the different plant tissues/pools (root, stem, leaf). DALEC-Grass can be integrated into a [model-data fusion (MDF) framework](https://www.sciencedirect.com/science/article/abs/pii/S0168192321001490) to constrain the dynamics of the C pools and fluxes using observed data. 

Stored in this repo are:
1. The original DALEC-Grass model code (DGv1)
2. The 2nd version of DG (DGv2) in which ACM2 has been integrated 
3. An experimental version of DG with N cycling 
4. Scipts to create DG inputs (DGIO), run the model data fusion (DGMDF) and forward run the DG (DGFWD) as well as a bash script to run each of these (cmd.sh)


![alt text](https://github.com/GCEL/DALEC-Grass/blob/master/dalec_grass.jpeg)
Figure: Schematic description of the DALEC-Grass model. DALEC-Grass simulates the dynamics of 5 C pools (C) : leaf, stem, roots, litter and SOC. C is allocated to the 5 C pools via NPP allocation (A) and litter production (L). Vegetation removals (VR) can occur due to grazing or cutting. DALEC-Grass determines whether a vegetation removal is caused by grazing, cutting or neither. When cutting is simulated (VRc>0) cut biomass (Bc) is removed from the ecosystem. When grazing is simulated (VRg>0) 32% of grazed biomass (Bg) is converted to manure, 54% of grazed biomass (Bg) is converted to animal respiration (aCO2) and 4%  of grazed biomass (Bg) is converted to methane (aCH4). Dotted lines show outward fluxes of C.  Solid lines show internal and inward fluxes of C



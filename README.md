# DALEC-GRASS
This fortran code is a development of the [Data Assimilation Linked Ecosystem Carbon (DALEC)](https://www.geos.ed.ac.uk/homes/mwilliam/DALEC.html) model that simulates carbon (C) dynamics in managed grassland ecosystems. DALEC requires 10 input variables and uses 33 parameters to estimate various C pools and fluxes. Model inputs, parameters, pools and fluxes are described in the first lines of the [code](https://github.com/vmyrgiotis/DALEC_Grass/blob/master/DALEC_GRASS.f90). The [Aggregated Canopy Model (ACM)](https://doi.org/10.1890/1051-0761(1997)007[0882:PGPPIT]2.0.CO;2) lies at the core of DALEC and is used to estimate daily Gross Primary Production (GPP). At each timestep, DALEC allocates the estimated Net Primary Production (NPP) to the different plant tissues/pools (root, stem, leaf). DALEC can be integrated into a [model-data fusion (MDF) framework](https://doi.org/10.1073/pnas.1515160113) to constrain the dynamics of the C pools and fluxes using observed data.

![alt text](https://github.com/GCEL/DALEC-Grass/blob/master/dalec_grass.gif)



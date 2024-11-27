#!/bin/bash

"""
[1] DGIO.py -> Input Output 
* Creates inputs for DALEC-GRASS (DG). 
* Requires spatial TS on meteorology (ERA5), LAI (Sentinel2) and SAR VV/VH backscatter (Sentinel1)
* Arguments : 'siteName' 'YYYYMMDD' 'YYYYMMDD' ParPostSample timeStep s1s1FusionBBox 

[2] DGMDF.py -> Model Data Fusion  
* Performs the Model Data Fusion (MDF) process.
* Requires the output of the DGIO.py -> SiteName_IO.nc
* Arguments :  'insFile' 'assimMetric' 'assimMeth' 'mdfMode'

[3] DGFWD.py -> Forward Runs  
* Performs forward runs wth DG. 
* Requires the parameter posteriors produced by DGMDF.py -> METHOD_METRIC_PostPars.csv 
* Arguments : 'insFile' 'parsPostFile' 

[*] Command line : bash cmd.sh prog method metric parallel
"""

### Set paths to code and data DIRECTORIES  
srcDir="Users/vm/Documents/Scripts/DALEC_Grass"
dataDir="Users/vm/Documents/Scripts/DALEC_Grass"

### Set the list of SITES 
sites=("GreatField" "Burrows" "DairyField")

### Set the python CODE that will be run 
prog="$1" # IO - MDF - FWD

### Set METRIC and METHOD for MDF 
method="$2" # DEMCZ - SA - MCMC (and DREAM - LHS - DREAM - SCEUA)
metric="$3" # ACCU - AREA - RMSE   

### Set PARALLELISATION of the MDF prog on/off
parallel="$4" # PRL for parallelised - SEQ for sequential implementation - NA for forward runs 

### Activate conda environment 
conda activate py3

### Compile Dalec-Grass .f90 
# f2py -c /${srcDir}/DGv2.f90 -m DG 

### Loop through all sites and run prog(s)
for field in "${sites[@]}"; do
	
	### Create DG inputs 
	if [[ "$prog" == "IO" ]]; then
			python /${srcDir}/DGIO.py "$field" "20170101" "20230101" 500 7 100
	fi
	
	### Run model-data fusion 
	if [[ "$prog" == "MDF" ]]; then
			if [[ "$parallel" == "PRL" ]]; then
				mpirun -n 6 python /${srcDir}/DGMDF.py "/${dataDir}/${field}/${field}_Inputs.nc" "${method}" "${metric}" "${parallel}"
			else
				python /${srcDir}/DGMDF.py "/${dataDir}/${field}/${field}_Inputs.nc" "${method}" "${metric}" "${parallel}"
			fi
	fi 
	
	### Run DG forward  
	if [[ "$prog" == "FWD" ]]; then
			python /${srcDir}/DGFWD.py "/${dataDir}/${field}/${field}_Inputs.nc" "/${dataDir}/${field}/${method}_${metric}_PostPars.csv"
	fi 

done
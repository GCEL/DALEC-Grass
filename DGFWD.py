# python SIM.py insFile parsPostFile 
import spotpy
import numpy as np 
import pandas as pd 
import DG
from sklearn.metrics import mean_squared_error
import xarray as xr
from math import sqrt
import math
from scipy import stats
import sys 
import netCDF4 as nc 
# import warnings
# warnings.filterwarnings('ignore')

### User input 
insFile = sys.argv[1]
parsPostFile = sys.argv[2]

### Posteriors Info 
prior_ranges = pd.read_csv(parsPostFile)
prior_ranges = prior_ranges.dropna()
prior_ranges = (prior_ranges.sort_values(by='like1',ascending=False))[:100]

### Drivers and Setup 
ins      = xr.open_dataset(insFile,engine='netcdf4')    
met      = np.vstack([ins.in_RunDay,ins.in_MinT,ins.in_MaxT,ins.in_sRad,ins.in_atmCO2,ins.in_DOY,ins.in_Prec,ins.in_VegMgmt,ins.in_BAF,ins.in_MinT21,ins.in_Photo21,ins.in_VPD21,ins.in_FMC,ins.in_avgT, ins.in_wind, ins.in_VPD])
deltat   = np.array(ins.in_deltat)
lat      = ins.attrs['simulation_site_lat']
nodays   = ins.attrs['simulation_timesteps_no'] 
start    = ins.attrs['simulation_start']   
finish   = ins.attrs['simulation_finish']   
spinup   = ins.attrs['simulation_spinup_days'] + 1  
nomet    = ins.attrs['simulation_met_no']    
nopools  = ins.attrs['simulation_pools_no']    
nofluxes = ins.attrs['simulation_fluxes_no']    
nopars   = ins.attrs['simulation_pars_no']    
nodaysOut = nodays - spinup # drop the 1st spinup year 

### RUN 
nosamples = len(prior_ranges)
unc_DF = pd.DataFrame()

### LAI GPP NEE
lai_mult = np.zeros([nosamples,nodaysOut])
gpp_mult = np.zeros([nosamples,nodaysOut])
nee_mult = np.zeros([nosamples,nodaysOut])

### FLUXES
autoResp_mult    = np.zeros([nosamples,nodaysOut]) # 2
respHetLit_mult  = np.zeros([nosamples,nodaysOut]) # 12
respHetSOM_mult  = np.zeros([nosamples,nodaysOut]) # 13
litter2SOM_mult  = np.zeros([nosamples,nodaysOut]) # 14

animmanure_mult = np.zeros([nosamples,nodaysOut]) # 18
animco2_mult = np.zeros([nosamples,nodaysOut])    # 19 
animch4_mult = np.zeros([nosamples,nodaysOut])    # 20
cut_mult = np.zeros([nosamples,nodaysOut])        # 21
graze_mult = np.zeros([nosamples,nodaysOut])      # 22 

evaporation_mult = np.zeros([nosamples,nodaysOut])   # 45 
transpiration_mult = np.zeros([nosamples,nodaysOut]) # 46
soilEvap_mult = np.zeros([nosamples,nodaysOut])      # 47 
wetCanEvap_mult = np.zeros([nosamples,nodaysOut])    # 48 
runoff_mult = np.zeros([nosamples,nodaysOut])        # 49
underflow_mult = np.zeros([nosamples,nodaysOut])     # 50 

### POOLS
stemC_mult = np.zeros([nosamples,nodaysOut+1])     # 0 
leafC_mult = np.zeros([nosamples,nodaysOut+1])     # 1
rootC_mult = np.zeros([nosamples,nodaysOut+1])     # 2
litterC_mult = np.zeros([nosamples,nodaysOut+1])   # 3
soilC_mult = np.zeros([nosamples,nodaysOut+1])     # 4 
soilWater_mult = np.zeros([nosamples,nodaysOut+1]) # 5


### Forward runs with N parameter vectors.
for i in range(0,nosamples):
	
	pars = np.array(prior_ranges[prior_ranges.columns[1:(nopars+1)]].iloc[i],order="F")
	
	lai,gpp,nee,fluxes,pools = DG.carbon_model_mod.carbon_model(start,finish,deltat,lat,met,pars,nopools,nofluxes)

	## Drop 1st/spinup year 
	lai = lai[spinup:]
	gpp = gpp[spinup:]
	nee = nee[spinup:]
	fluxes = fluxes[spinup:,:]
	pools = pools[spinup:,:]

	### LAI GPP NEE 
	lai_mult[i,:] = lai
	gpp_mult[i,:] = gpp
	nee_mult[i,:] = nee 
 
	### FLUXES
	autoResp_mult[i,:] =  fluxes[:,2]
	respHetLit_mult[i,:] =  fluxes[:,12]
	respHetSOM_mult[i,:] =  fluxes[:,13]
	litter2SOM_mult[i,:] =  fluxes[:,14]
	animmanure_mult[i,:] =  fluxes[:,18]
	animco2_mult[i,:] =  fluxes[:,19]
	animch4_mult[i,:] =  fluxes[:,20]
	cut_mult[i,:] =  fluxes[:,21]
	graze_mult[i,:] =  fluxes[:,22]
	evaporation_mult[i,:] =  fluxes[:,45]
	transpiration_mult[i,:] =  fluxes[:,46]
	soilEvap_mult[i,:] =  fluxes[:,47]
	wetCanEvap_mult[i,:] =  fluxes[:,48]
	runoff_mult[i,:] =  fluxes[:,49]
	underflow_mult[i,:] =  fluxes[:,50]

	### POOLS
	stemC_mult[i,:] = pools[:,0]
	leafC_mult[i,:] = pools[:,1]
	rootC_mult[i,:] = pools[:,2]
	litterC_mult[i,:] = pools[:,3]
	soilC_mult[i,:] = pools[:,4]
	soilWater_mult[i,:] = pools[:,5]

###  Estimate mean/SD for considered variables and save to unc_DF
for i in range(nodaysOut): 
	unc_DF = unc_DF.append({
			  ### LAI GPP NEE        
			  'lai': stats.bayes_mvs(lai_mult[:,i],alpha=0.99)[0][0] , 
			  'lai_std':stats.bayes_mvs(lai_mult[:,i],alpha=0.99)[2][0],
			  'lai_normality_metric': np.mean(lai_mult[:,i])/ np.median(lai_mult[:,i]) ,              
			  'gpp': stats.bayes_mvs(gpp_mult[:,i],alpha=0.99)[0][0] , 
			  'gpp_std': stats.bayes_mvs(gpp_mult[:,i],alpha=0.99)[2][0],
			  'gpp_normality_metric':  np.mean(gpp_mult[:,i])/ np.median(gpp_mult[:,i]) ,             
			  'nee': stats.bayes_mvs(nee_mult[:,i],alpha=0.99)[0][0] , 
			  'nee_std': stats.bayes_mvs(nee_mult[:,i],alpha=0.99)[2][0],
			  'nee_normality_metric':  np.mean(nee_mult[:,i])/ np.median(nee_mult[:,i]) ,                
			  ### FLUXES
			  'autoResp': stats.bayes_mvs(autoResp_mult[:,i],alpha=0.99)[0][0] , 
			  'autoResp_std':  stats.bayes_mvs(autoResp_mult[:,i],alpha=0.99)[2][0],    
			  'respHetLit': stats.bayes_mvs(respHetLit_mult[:,i],alpha=0.99)[0][0] , 
			  'respHetLit_std':  stats.bayes_mvs(respHetLit_mult[:,i],alpha=0.99)[2][0],    
			  'respHetSOM': stats.bayes_mvs(respHetSOM_mult[:,i],alpha=0.99)[0][0] , 
			  'respHetSOM_std':  stats.bayes_mvs(respHetSOM_mult[:,i],alpha=0.99)[2][0],    
			  'litter2SOM': stats.bayes_mvs(litter2SOM_mult[:,i],alpha=0.99)[0][0] , 
			  'litter2SOM_std':  stats.bayes_mvs(litter2SOM_mult[:,i],alpha=0.99)[2][0],    
			  'animManure': stats.bayes_mvs(animmanure_mult[:,i],alpha=0.99)[0][0] , 
			  'animManure_std':  stats.bayes_mvs(animmanure_mult[:,i],alpha=0.99)[2][0],    
			  'animCO2': stats.bayes_mvs(animco2_mult[:,i],alpha=0.99)[0][0] , 
			  'animCO2_std':  stats.bayes_mvs(animco2_mult[:,i],alpha=0.99)[2][0],    
			  'animCH4': stats.bayes_mvs(animch4_mult[:,i],alpha=0.99)[0][0] , 
			  'animCH4_std':  stats.bayes_mvs(animch4_mult[:,i],alpha=0.99)[2][0],    
			  'cut': stats.bayes_mvs(cut_mult[:,i],alpha=0.99)[0][0] , 
			  'cut_std':  stats.bayes_mvs(cut_mult[:,i],alpha=0.99)[2][0],    
			  'graze': stats.bayes_mvs(graze_mult[:,i],alpha=0.99)[0][0] , 
			  'graze_std':  stats.bayes_mvs(graze_mult[:,i],alpha=0.99)[2][0],    
			  'evaporation': stats.bayes_mvs(evaporation_mult[:,i],alpha=0.99)[0][0] , 
			  'evaporation_std':  stats.bayes_mvs(evaporation_mult[:,i],alpha=0.99)[2][0],    
			  'transpiration': stats.bayes_mvs(transpiration_mult[:,i],alpha=0.99)[0][0] , 
			  'transpiration_std':  stats.bayes_mvs(transpiration_mult[:,i],alpha=0.99)[2][0],    
			  'soilEvap': stats.bayes_mvs(soilEvap_mult[:,i],alpha=0.99)[0][0] , 
			  'soilEvap_std':  stats.bayes_mvs(soilEvap_mult[:,i],alpha=0.99)[2][0],    
			  'wetCanEvap': stats.bayes_mvs(wetCanEvap_mult[:,i],alpha=0.99)[0][0] , 
			  'wetCanEvap_std':  stats.bayes_mvs(wetCanEvap_mult[:,i],alpha=0.99)[2][0],    
			  'runoff': stats.bayes_mvs(runoff_mult[:,i],alpha=0.99)[0][0] , 
			  'runoff_std':  stats.bayes_mvs(runoff_mult[:,i],alpha=0.99)[2][0],    
			  'underflow': stats.bayes_mvs(underflow_mult[:,i],alpha=0.99)[0][0] , 
			  'underflow_std':  stats.bayes_mvs(underflow_mult[:,i],alpha=0.99)[2][0],    
			  ### POOLS
			  'leafC': stats.bayes_mvs(leafC_mult[:,1+i],alpha=0.99)[0][0] , 
			  'leafC_std': stats.bayes_mvs(leafC_mult[:,1+i],alpha=0.99)[2][0],
			  'stemC': stats.bayes_mvs(stemC_mult[:,1+i],alpha=0.99)[0][0] , 
			  'stemC_std': stats.bayes_mvs(stemC_mult[:,1+i],alpha=0.99)[2][0],
			  'rootC': stats.bayes_mvs(rootC_mult[:,1+i],alpha=0.99)[0][0] , 
			  'rootC_std': stats.bayes_mvs(rootC_mult[:,1+i],alpha=0.99)[2][0],
			  'litterC': stats.bayes_mvs(litterC_mult[:,1+i],alpha=0.99)[0][0] , 
			  'litterC_std': stats.bayes_mvs(litterC_mult[:,1+i],alpha=0.99)[2][0],  
			  'soilC': stats.bayes_mvs(soilC_mult[:,1+i],alpha=0.99)[0][0] , 
			  'soilC_std': stats.bayes_mvs(soilC_mult[:,1+i],alpha=0.99)[2][0],
			  'SW': stats.bayes_mvs(soilWater_mult[:,1+i],alpha=0.99)[0][0] , 
			  'SW_std': stats.bayes_mvs(soilWater_mult[:,1+i],alpha=0.99)[2][0],
	},ignore_index=True)


## Add predicted LAI GPP NEE Pools Fluxes and save to netcdf
### LAI GPP NEE
ins.LAI.data=np.array(unc_DF.lai)
ins.LAI_sd.data=np.array(unc_DF.lai_std)
ins.GPP.data=np.array(unc_DF.gpp)
ins.GPP_sd.data=np.array(unc_DF.gpp_std)
ins.NEE.data=np.array(unc_DF.nee)
ins.NEE_sd.data=np.array(unc_DF.nee_std)
### FLUXES
ins.flux_autoResp.data=np.array(unc_DF.autoResp)
ins.flux_autoResp_sd.data=np.array(unc_DF.autoResp_std)
ins.flux_hetResLit.data=np.array(unc_DF.respHetLit)
ins.flux_hetResLit_sd.data=np.array(unc_DF.respHetLit_std)
ins.flux_hetResSOM.data=np.array(unc_DF.respHetSOM)
ins.flux_hetResSOM_sd.data=np.array(unc_DF.respHetSOM_std)
ins.flux_lit2SOM.data=np.array(unc_DF.litter2SOM)
ins.flux_lit2SOM_sd.data=np.array(unc_DF.litter2SOM_std)
ins.flux_animal_manure_to_soil.data=np.array(unc_DF.animManure)
ins.flux_animal_manure_to_soil_sd.data=np.array(unc_DF.animManure_std)
ins.flux_animal_respiration.data=np.array(unc_DF.animCO2)
ins.flux_animal_respiration_sd.data=np.array(unc_DF.animCO2_std)
ins.flux_animal_methane.data=np.array(unc_DF.animCH4)
ins.flux_animal_methane_sd.data=np.array(unc_DF.animCH4_std)
ins.flux_harvest.data=np.array(unc_DF.cut)
ins.flux_harvest_sd.data=np.array(unc_DF.cut_std)
ins.flux_grazed.data=np.array(unc_DF.graze)
ins.flux_grazed_sd.data=np.array(unc_DF.graze_std)
ins.flux_evap.data=np.array(unc_DF.evaporation)
ins.flux_evap_sd.data=np.array(unc_DF.evaporation_std)
ins.flux_transpir.data=np.array(unc_DF.transpiration)
ins.flux_transpir_sd.data=np.array(unc_DF.transpiration_std)
ins.flux_soilEvap.data=np.array(unc_DF.soilEvap)
ins.flux_soilEvap_sd.data=np.array(unc_DF.soilEvap_std)
ins.flux_wetCanopEvap.data=np.array(unc_DF.wetCanEvap)
ins.flux_wetCanopEvap_sd.data=np.array(unc_DF.wetCanEvap_std)
ins.flux_runoff.data=np.array(unc_DF.runoff)
ins.flux_runoff_sd.data=np.array(unc_DF.runoff_std)
ins.flux_underflow.data=np.array(unc_DF.underflow)
ins.flux_underflow_sd.data=np.array(unc_DF.underflow_std)
### POOLS
ins.pool_leaf.data=np.array(unc_DF.leafC)
ins.pool_leaf_sd.data=np.array(unc_DF.leafC_std)
ins.pool_stem.data=np.array(unc_DF.stemC)
ins.pool_stem_sd.data=np.array(unc_DF.stemC_std)
ins.pool_root.data=np.array(unc_DF.rootC)
ins.pool_root_sd.data=np.array(unc_DF.rootC_std)
ins.pool_SOM.data=np.array(unc_DF.soilC)
ins.pool_SOM_sd.data=np.array(unc_DF.soilC_std)
ins.pool_litter.data=np.array(unc_DF.litterC)
ins.pool_litter_sd.data=np.array(unc_DF.litterC_std)
ins.pool_SW.data=np.array(unc_DF.SW)
ins.pool_SW_sd.data=np.array(unc_DF.SW_std)
### SAVE NETCDF 
ins.to_netcdf(insFile)
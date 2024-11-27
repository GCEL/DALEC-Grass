"""
- python MDF.py insFile assimMetric assimMeth mdfMode
- always check against DGIO code and outputs 
"""
import pandas as pd 
import numpy as np
import spotpy
from sklearn.metrics import mean_squared_error
from math import sqrt
import DG
import xarray as xr 
import sys 
from scipy.integrate import simpson
import warnings
warnings.filterwarnings("ignore")

""" USER INPUT """
insFile = str(sys.argv[1]) # path to input data file 
assimMeth = str(sys.argv[2]) # method to use for assimilating obs data 
assimMetric = str(sys.argv[3]) # metric to use for obs data assimilation  
mdfMode = str(sys.argv[4]) # run assimilation method using parallelisation or not (on/off)
workDir = '/'
for i in range(len(insFile.split('/'))-1) : workDir = workDir + str(insFile.split('/')[i]) + '/'

class abc_dalec() :

	def __init__(self):   
		
		""" DRIVERS AND SETUP INFO """ 
		self.ins      = xr.open_dataset(insFile,engine='netcdf4')  
		self.met      = np.vstack([self.ins.in_RunDay,self.ins.in_MinT,self.ins.in_MaxT,self.ins.in_sRad,self.ins.in_atmCO2,self.ins.in_DOY,self.ins.in_Prec,self.ins.in_VegMgmt,
								   self.ins.in_BAF,self.ins.in_MinT21,self.ins.in_Photo21,self.ins.in_VPD21,self.ins.in_FMC,self.ins.in_avgT, self.ins.in_wind, self.ins.in_VPD])
		self.lat      = self.ins.attrs['simulation_site_lat']
		self.deltat   = np.array(self.ins.in_deltat)
		self.nodays   = int(self.ins.attrs['simulation_timesteps_no'])
		self.start    = int(self.ins.attrs['simulation_start'])
		self.finish   = int(self.ins.attrs['simulation_finish'])
		self.spinup   = int(self.ins.attrs['simulation_spinup_days']) + 1 
		self.nomet    = int(self.ins.attrs['simulation_met_no'])
		self.nopools  = int(self.ins.attrs['simulation_pools_no'])
		self.nofluxes = int(self.ins.attrs['simulation_fluxes_no'])
		self.nopars   = int(self.ins.attrs['simulation_pars_no'])
		self.met[7,:self.spinup] = 0 ### Free growing grass during spinup

		### OBSERVATIONAL DATA
		self.laiDF           = pd.DataFrame(data=self.ins.obs_LAI[self.spinup:], index=pd.to_datetime(np.array(self.ins.obs_LAI_dates[self.spinup:])), columns=['obs'])
		self.laiDF['obsUnc'] = 2*self.ins.obs_LAI_sd[self.spinup:] # obs LAI uncertainty (relative SD)
		self.cutDF           = pd.DataFrame(index=pd.to_datetime(np.array(self.ins.obs_LAI_dates[self.spinup:])), columns=['sim'])

		### OBSERVATIONAL DATA
		self.neeDF           = pd.DataFrame(data=self.ins.obs_NEE[self.spinup:], index=pd.to_datetime(np.array(self.ins.obs_NEE_dates[self.spinup:])), columns=['obs'])
        
		# # ### REFINED PARAMETER PRIORS 
		# self.pars = pd.read_csv("/data/notebooks/jupyterlab-0001/DALEC_GRASS/data/SiteData/GreatField/RefPriors.csv")    
		# self.pars = self.pars.sort_values(by='like1',ascending=False)[:100]

	
	""" MODEL PARAMETER PRIOR RANGES """ 
	def parameters(self): 
	
		self.params=[]

		for x in range(0,self.nopars) : 
			### INPUT FILE PRIORS 
			self.params.append(spotpy.parameter.Uniform(str(f'P{x+1}'), self.ins[f'par_P{x+1}'].prior_min, self.ins[f'par_P{x+1}'].prior_max)) # edited priors 
			# self.params.append(spotpy.parameter.Uniform(str(f'P{x+1}'), self.ins[f'par_P{x+1}'].def_prior_min, self.ins[f'par_P{x+1}'].def_prior_max)) # def priors - Luke's
			### OTHER REFINED PRIORS 
			# self.params.append(spotpy.parameter.Uniform('P%s' %(x+1), self.pars['parP%s'%(x+1)].min(), self.pars['parP%s'%(x+1)].max())) # previously produced/refined priors
																						
		return spotpy.parameter.generate(self.params)

	
	""" MODEL IMPLEMENTATION """
	def simulation(self,vector):   

		pars = np.array(vector, order='F')
		
		### PARAMETER SANITY CHECK 
		if ( (pars[26]>pars[27]) & (pars[12]<=pars[11]) or (pars[21]<=pars[20]) or (pars[19]<=pars[13]) or (pars[22]<pars[15:19].sum()) or 
            (pars[17]>pars[15]+pars[16]+pars[18]) or (pars[17]/(pars[15]+pars[16])<3) or (pars[7]>pars[6]) ) : return [-np.inf]

		else :
						
			### Run DALEC-GRASS
			lai,gpp,nee,fluxes,pools = DG.carbon_model_mod.carbon_model(self.start,self.finish,self.deltat,self.lat,self.met,pars,self.nopools,self.nofluxes)
					
			### Drop spinup year from outputs 
			lai = lai[self.spinup:]
			gpp = gpp[self.spinup:]
			nee = nee[self.spinup:]
			fluxes = fluxes[self.spinup:,:]
			pools = pools[self.spinup:,:]

			### Add sim LAI to obs LAI DF 
			self.laiDF['sim'] = lai
			self.laiDF['relError'] = abs(self.laiDF.obs - self.laiDF.sim) / self.laiDF.obs # difference between sim and obs LAI (relative to obs)
			sub_laiDF = self.laiDF[self.laiDF.obs >= 1] # Drop LAI obs < 1 m2.m-2

			### Add sim NEE to obs LAI NEE 
			self.neeDF['sim'] = nee
			sub_neeDF = self.neeDF.dropna()
            
			### Calculate annual LAI peaks 
			annual_laiDF_peaks = sub_laiDF.resample('YE').max() 
            
			### DF to hold annual sum info on grass cutting / harvest 
			self.cutDF['sim'] = fluxes[:,21]
			sub_cutDF = self.cutDF.sim.resample('YE').sum()
			
			### ECOLOGICAL AND DYNAMIC CONSTRAINTS 
			if  ((np.isnan(pools).any()) or (np.isnan(lai).any()) or (np.isnan(gpp).any()) or 
				(np.any(pools < 0)) or (np.any(fluxes < 0)) or (np.any(lai < 0)) or (np.all(pools == 0, axis=0)).any() or     
				(np.any(fluxes[:,21]*self.deltat[0] * 0.021 > 10)) or # havest yield cannot be> X t.DM.ha-1)
				(np.any(fluxes[:,22]*self.deltat[0] * 0.021 > 10)) or # grazed biomass in a week/timestep cannot be > X tDMha-1
				((sub_cutDF==0).all()) or  # at least one cut over the simulated period 
				(abs((pools[-1,4]/pools[0,4])-1) > 0.05) or # ~ stable SOC pool
				(np.any(annual_laiDF_peaks.sim/annual_laiDF_peaks.obs < 0.8)) # capture annual LAI peaks (mean lAI RelUnc: 20%)
				# (round(sqrt(mean_squared_error(sub_laiDF.obs,sub_laiDF.sim)),2) > 1) # obs-vs-sim LAI RMSE > 1 
				) : return [-np.inf]
				
		""" CHOOOSE THE SIM-VS-OBS METRIC """
		if assimMetric == 'AREA' : # AREA under LAI curve (0-inf)
			return [abs(1-round(simpson(np.array(sub_laiDF.sim),dx=5)/simpson(np.array(sub_laiDF.obs),dx=5),2))]
		
		if assimMetric == 'ACCU' : # ACCUracy metric (0-1) 
			return [round((len(sub_laiDF[sub_laiDF.relError <= sub_laiDF.obsUnc])/len(sub_laiDF)),2)] 

		if assimMetric == 'ARAC' : # AREA AND ACCU 
			return [[abs(1-round(simpson(np.array(sub_laiDF.sim),dx=5)/simpson(np.array(sub_laiDF.obs),dx=5),3)),
					round((len(sub_laiDF[sub_laiDF.relError <= sub_laiDF.obsUnc])/len(sub_laiDF)),2)]]

		if assimMetric == 'ARRM' : # AREA AND ACCU 
			return [np.mean([abs(1-round(simpson(np.array(sub_laiDF.sim),dx=5)/simpson(np.array(sub_laiDF.obs),dx=5),3)),
							round(sqrt(mean_squared_error(sub_laiDF.obs,sub_laiDF.sim)),2)])]
        
		if assimMetric == 'RMSE' : # Root Mean Squared Error 
			# return [round(sqrt(mean_squared_error(sub_laiDF.obs,sub_laiDF.sim)),2)] 
			return [round(sqrt(mean_squared_error(sub_neeDF.obs,sub_neeDF.sim)),2)] 

	""" OBJECTIVE FUNCTION TARGET """
	def evaluation(self) : 
		
		if assimMetric == 'AREA' : return [0]     # when using AREA -> 0 means sim and obs area-under-curve are equal
		if assimMetric == 'ACCU' : return [1]     # when using accuracy -> 1 means perfect obs-vs-sim fit 		
		if assimMetric == 'ARAC' : return [[0,1]] # when using AREA + ACCU
		if assimMetric == 'ARRM' : return [0]     # when using AREA + RMSE
		if assimMetric == 'RMSE' : return [0]     # when using RMSE -> 0 means perfect obs-vs-sim fit 

	""" OBJECTIVE FUNCTION FORMULA """
	def objectivefunction(self,simulation,evaluation) :
		objectivefunction = -spotpy.objectivefunctions.mae(evaluation,simulation)    
		return objectivefunction


""" CHOOSE ASSIMILATION/SAMPLING ALGORITHM AND PARALLILISATION """
if mdfMode == 'PRL' : 
	
	### MCMC Metropolis-Hastings (MCMC)
	if assimMeth == 'MCMC':
		results = [] 
		spotpy_setup = abc_dalec()     
		sampler = spotpy.algorithms.mcmc(spotpy_setup, dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars', dbformat='csv', save_sim=True, parallel='mpi')
		results.append(sampler.sample(10000000,nChains=10))

	### Differential Evolution Markov Chain (DEMCZ)
	if assimMeth == 'DEMCZ':
		results = [] 
		spotpy_setup = abc_dalec()     
		sampler = spotpy.algorithms.demcz(spotpy_setup, dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars', dbformat='csv', save_sim=True, parallel='mpi')
		results.append(sampler.sample(100000, nChains=10))

	### Simulated Annealing (SA)
	if assimMeth == 'SA':
		results = [] 
		spotpy_setup = abc_dalec() 
		sampler = spotpy.algorithms.sa(spotpy_setup,dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars',dbformat='csv',  save_sim=True, parallel='mpi')
		results.append(sampler.sample(repetitions=100000000,Tini=98,Ntemp=1000,alpha=0.99)) # tini: Starting temperature | Ntemp: No of trials per T | alpha: T reduction


if mdfMode == 'SEQ' : 

	### Latin Hypercube Sampling (LHS)
	if assimMeth == 'LHS':
		results = [] 
		spotpy_setup = abc_dalec()         
		sampler = spotpy.algorithms.lhs(spotpy_setup, dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars',  dbformat='csv', save_sim=True)
		results.append(sampler.sample(repetitions=10000000))

	### Differential Evolution Adaptive Metropolis Algorithm (DREAM)
	if assimMeth == 'DREAM':
		results = [] 
		spotpy_setup = abc_dalec()         
		sampler = spotpy.algorithms.dream(spotpy_setup, dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars',  dbformat='csv', save_sim=True)
		results.append(sampler.sample(repetitions=1000000,nChains=3*37+1))

	### Shuffled Complex Evolution Algorithm (SCEUA)
	if assimMeth == 'SCEUA':
		results = [] 
		spotpy_setup = abc_dalec()     
		sampler = spotpy.algorithms.sceua(spotpy_setup, dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars', dbformat='csv', save_sim=True)
		results.append(sampler.sample(repetitions=10000000,ngs=2*37))

	### MCMC Metropolis-Hastings (MCMC)
	if assimMeth == 'MCMC':
		results = [] 
		spotpy_setup = abc_dalec()     
		sampler = spotpy.algorithms.mcmc(spotpy_setup, dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars', dbformat='csv', save_sim=True)
		results.append(sampler.sample(10000000,nChains=10))

	### Differential Evolution Markov Chain (DEMCZ)
	if assimMeth == 'DEMCZ':
		results = [] 
		spotpy_setup = abc_dalec()     
		sampler = spotpy.algorithms.demcz(spotpy_setup, dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars', dbformat='csv', save_sim=True)
		results.append(sampler.sample(100000, nChains=100))

	### Simulated Annealing (SA)
	if assimMeth == 'SA':
		results = [] 
		spotpy_setup = abc_dalec() 
		sampler = spotpy.algorithms.sa(spotpy_setup,dbname=f'{workDir}{assimMeth}_{assimMetric}_PostPars',dbformat='csv', save_sim=True)
		results.append(sampler.sample(repetitions=100000000,Tini=98,Ntemp=1000,alpha=0.99)) # tini: Starting temperature | Ntemp: No of trials per T | alpha: T reduction        
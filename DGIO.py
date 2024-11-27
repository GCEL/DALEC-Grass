"""
* Create a netcdf that contains model inputs, metadata and model outputs 
* command line : python DGv2_IO.py 'site' firstDate lastDate posteriorLen timestep s1s1fusionbbox 

* python DGv2_IO.py 'Burrows' 20170101 20230101 1000 7 100 
* 'GreatField' : 2 'Burrows' : 4 'DairyField' :  9 

! Always check met data source ERA5 vs ground measured 
"""
import os ; os.chdir('/Users/vm/Documents/Scripts/DALEC_Grass/')
from auxFunc import daylength as DC
from auxFunc import S1S2_Fusion as LAICalc
from auxFunc import SoilGridsSOC as SOCCalc
from shapely import wkt
from datetime import datetime, timedelta
import json
import netCDF4 as nc
from glob import glob 
import pandas as pd
import geopandas as gpd 
import numpy as np 
import xarray as xr 
import subprocess
import sys 
import math
import warnings
# warnings.filterwarnings("ignore")


"""
INFO ON DIRECTORIES AND CODE INPUTS 
"""

workdir = '/Users/vm/Desktop/dir'

# ## Input info or testing 
# site='GreatField'
# sitePoly = gpd.read_file(f'{workdir}/RRNW/{site}.geojson') 
# firstDate="20170101" 
# lastDate="20230101" 
# posteriorLen=500 
# timestep=7 
# S1S2BBoxSize=50

### Command line / program user inputs 
site = str(sys.argv[1]) # str
sitePoly = gpd.read_file(f'{workdir}/RRNW/{site}.geojson') 
ClipPoly = sitePoly.geometry.to_list() 
firstDate = str(sys.argv[2]) # str YYYYMMDD
lastDate = str(sys.argv[3])  # str YYYYMMDD
posteriorLen = int(sys.argv[4]) # int 
timestep = int(sys.argv[5]) # int (1 or 7) 
S1S2BBoxSize = int(sys.argv[6]) # int 
# S1FilesPath = str(sys.argv[7]) # str /path/to/dir 
# S2FilesPath = str(sys.argv[8]) # str /path/to/dir
# ERA5FilesPath = str(sys.argv[9]) # str /path/to/dir



"""
SOIL C POOL SIZE ESTIMATE
"""

# soilCValue = round(SOCCalc(f'{workdir}/RRNW/SoilGrids/',f'{workdir}/RRNW/{site}.geojson'))

if site == 'GreatField' : soilCValue = 1776
if site == 'DairyField' : soilCValue = 1774
if site == 'Burrows'    : soilCValue = 2734

"""
CREATE MET DRIVERS TIME SERIES 
"""

### GROUND MEASURED MET DATA 
##############################

# ### Produce met_DF dataset 
# folders_met = []
# for file in glob("/Users/vm/Desktop/dir.nosync/RRNW/ground_data/measurements*") : 
# 	folders_met.append(file)
# folders_met.sort()

# df = pd.DataFrame()
# for F in folders_met : df = df._append(pd.read_csv(F))

# df.index = pd.to_datetime(df.Datetime)
# df = df[df.columns[df.columns.str.contains('Site')]]
# df = df[df.columns[~df.columns.str.contains('Quality')]]

# ## vapor pressure deficit (http:/cronklab.wikidot.com/calculation-of-vapour-pressure-deficit)
# df['VPD'] = (1-(df['Relative Humidity (% RH) [Site]']/100)) * (610.7*10**(7.5*df['Air Temperature (oC) [Site]']/(237.3+df['Air Temperature (oC) [Site]'])))

# met_DF = pd.DataFrame(columns=['date','minT','maxT','srad','vpd','wind','photopd','21d_vpd','21d_minT','21d_photopd','prec'])

# siteLat = round(gpd.read_file(f'{workdir}/RRNW/{site}.geojson').geometry.centroid.y[0],2)

# for t in range(len(df.resample('D').mean())):
# 	met_DF = met_DF._append({'date' : (df.resample('D').mean()).index[t],
# 							'doy'  : (df.resample('D').mean()).index.dayofyear[t],
# 							'minT' : (df['Air Temperature (oC) [Site]'].resample('D',label='left').min()).iloc[t],
# 							'maxT' : (df['Air Temperature (oC) [Site]'].resample('D',label='left').max()).iloc[t],
# 							'avgT' : (df['Air Temperature (oC) [Site]'].resample('D',label='left').mean()).iloc[t],
# 							'prec' : (df['Precipitation (mm) [Site]'].resample('D',label='left').mean()).iloc[t] / 900, # mm to mm.s-1 (i.e. kg.rain.m-2.s-1) 900secs in 15mins 
# 							'srad' : (df['Solar Radiation (W/m2) [Site]'].resample('D',label='left').sum()).iloc[t] * 0.0036 / 4 , # raw data come very 15min W/m2 (per hour) * 0.0036 --> MJ.m-2.d-1
# 							'wind' : (df['Wind Speed (km/h) [Site]'].resample('D',label='left').mean()).iloc[t] * 0.277778, # km/h to m/s
# 							'vpd'  : (df.VPD.resample('D',label='left').mean()).iloc[t],
# 							'photopd':  DC((df.resample('D').mean()).index.dayofyear[t],float(siteLat))
# 							},ignore_index=True)

# ## 21-day rolling average photoperiod - minT - VPD
# met_DF['21d_vpd'] = met_DF['vpd'].rolling(window=21).mean() 
# met_DF['21d_vpd'] = met_DF['21d_vpd'].fillna(method='bfill')
# met_DF['21d_minT'] = met_DF['minT'].rolling(window=21).mean() + 273.15 # in Kelvin ?
# met_DF['21d_minT'] = met_DF['21d_minT'].fillna(method='bfill')
# met_DF['21d_photopd'] = met_DF['photopd'].rolling(window=21).mean() * 3600 # hrs to sec 
# met_DF['21d_photopd'] = met_DF['21d_photopd'].fillna(method='bfill')
# met_DF.to_csv(f'{workdir}/RRNW/RRNW_MetData.csv')


### ERA5 MET DATA 
###########################

## Download data (if not there)
if len(glob(f"{workdir}/ERA5/*.nc")) == 0 :
	for yyyy in np.arange(2015,2023) : 
		auxFunc.ERA5_DL(year=yyyy, path2outs=f"{workdir}/ERA5/",
						AOIPolygon=[sitePoly.geometry.centroid.y[0]+0.1,
									sitePoly.geometry.centroid.x[0]-0.1,
									sitePoly.geometry.centroid.y[0]-0.1,
									sitePoly.geometry.centroid.x[0]+0.1])

## Produce met_DF dataset 
folders_met = []
for file in glob(f"{workdir}/ERA5/*.nc") : 
	if int(file.split('_')[1][0:4]) in np.arange(int(firstDate[:4])-1,int(lastDate[:4])) : # CHECK , add 1yr of spinup met inputs 
		folders_met.append(file)
folders_met.sort()

siteLat = sitePoly.geometry.centroid.y[0] #.centroid.y.iloc[0]
siteLon = sitePoly.geometry.centroid.x[0] #.centroid.x.iloc[0]

met_DF = pd.DataFrame(columns=['date','minT','maxT','srad','vpd','wind','photopd','21d_vpd','21d_minT','21d_photopd','precipitation'])

for ii in range(len(folders_met)):

	ds = xr.open_dataset(folders_met[ii])
	dsloc = ds.sel(longitude=float(siteLon),latitude=float(siteLat),method='nearest') 
	df = dsloc.to_dataframe()
	## unit conversion
	df.t2m = df.t2m - 273.15
	df.d2m = df.d2m - 273.15
	## relative humidity (http:/andrew.rsmas.miami.edu/bmcnoldy/Humidity.html)
	df['RH'] = 100*(np.exp((17.625*df.d2m)/(243.04+df.d2m))/np.exp((17.625*df.t2m)/(243.04+df.t2m))) 
	## vapor pressure deficit (http:/cronklab.wikidot.com/calculation-of-vapour-pressure-deficit)
	df['VPD'] = (1-(df.RH/100)) * (610.7*10**(7.5*df.t2m/(237.3+df.t2m)))
	df['windspeed'] = (df.u10**2 + df.v10**2)**0.5

	df = df.reset_index()
	df.index = df.time

	for t in range(len(df.resample('D').mean())):
			met_DF = met_DF._append({'date' : (df.resample('D').mean()).index[t],
									'doy'  : (df.resample('D').mean()).index.dayofyear[t],
									'minT' : (df.t2m.resample('D',label='left').min()).iloc[t],
									'maxT' : (df.t2m.resample('D',label='left').max()).iloc[t],
									'avgT' : (df.t2m.resample('D',label='left').mean()).iloc[t],
									'prec' : (df.tp.resample('D',label='left').mean()).iloc[t] * 1000 / 3600 , # m to mm.s-1 i.e. kg.rain.m-2.s-1
									'srad' : (df.ssrd.resample('D',label='left').sum()).iloc[t] * 1e-6, # to MJ.m-2.d-1
									'wind' : (df.windspeed.resample('D',label='left').mean()).iloc[t],
									'vpd'  : (df.VPD.resample('D',label='left').mean()).iloc[t],
									'photopd':  DC((df.resample('D').mean()).index.dayofyear[t],float(siteLat))
									},ignore_index=True)

	## 21-day rolling average photoperiod - minT - VPD
	met_DF['21d_vpd'] = met_DF['vpd'].rolling(window=21).mean() 
	met_DF['21d_vpd'] = met_DF['21d_vpd'].fillna(method='bfill')
	met_DF['21d_minT'] = met_DF['minT'].rolling(window=21).mean() + 273.15
	met_DF['21d_minT'] = met_DF['21d_minT'].fillna(method='bfill')
	met_DF['21d_photopd'] = met_DF['photopd'].rolling(window=21).mean() * 3600 # hrs to sec 
	met_DF['21d_photopd'] = met_DF['21d_photopd'].fillna(method='bfill')

### Produce DG inputs from met_DF
# met_DF = pd.read_csv(f'{workdir}/RRNW/RRNW_MetData.csv') # if already produced and saved as .csv
met_DF.index = pd.to_datetime(met_DF.date)
met_DF = met_DF[(met_DF.index>=str(int(firstDate)-10000))&(met_DF.index<=lastDate)] # load met_DF and keep data for period : (firstDate - spinupYear) - (lastDate)
met_DF = met_DF.drop(columns=['date','photopd'])

# ### T, SRAD, VPD, Photoperiod time-series # CHECK 
weekly_tmax = met_DF.maxT.resample('7D',label='left').mean()
weekly_tmin = met_DF.minT.resample('7D',label='left').mean()
weekly_tavg = met_DF.avgT.resample('7D',label='left').mean()
weekly_prec = met_DF.prec.resample('7D',label='left').mean()
weekly_rad  = met_DF['srad'].resample('7D',label='left').mean()
weekly_DOY  = met_DF.doy.resample('7D',label='left').max()
weekly_21d_vpd = met_DF['21d_vpd'].resample('7D',label='left').mean()
weekly_21d_minT = met_DF['21d_minT'].resample('7D',label='left').mean()
weekly_21d_photopd = met_DF['21d_photopd'].resample('7D',label='left').mean()
weekly_vpd  = met_DF['vpd'].resample('7D',label='left').mean()
weekly_wind  = met_DF['wind'].resample('7D',label='left').mean()

weekly_in_dates = np.array([])
for i in range(len(weekly_21d_photopd)) : 
	weekly_in_dates = np.append(str(weekly_21d_photopd.index.date[i]),weekly_in_dates) 
weekly_in_dates = np.flipud(weekly_in_dates)

### Atmospheric CO2 TS - global mean from Mauna Loa
subprocess.call(f"rm {workdir}/weekly_in_situ_co2_mlo.csv",shell=True)
subprocess.call(f"wget https://scrippsco2.ucsd.edu/assets/data/atmospheric/stations/in_situ_co2/weekly/weekly_in_situ_co2_mlo.csv -P {workdir}",shell=True)
co2 = pd.read_csv(f'{workdir}/weekly_in_situ_co2_mlo.csv', skiprows=44, na_values=["", "NA"], header=None, keep_default_na=False)
co2 = co2.rename(columns={0:'date',1:'ppm'})
co2["date"] = pd.to_datetime(co2["date"], format="%Y-%m-%d")
co2 = co2.set_index("date")
co2_ppm = co2.ppm.resample('7D',label='left').mean() # CHECK 
co2_ppm[co2_ppm<=0] = np.nan 
co2_ppm = co2_ppm.interpolate()
co2_ppm = co2_ppm[pd.to_datetime(co2_ppm.index.date,format='%Y-%m-%d') >= weekly_21d_photopd.index[0]][:(len(weekly_21d_minT))] # CHECK , ensure spinup year data are added


"""
CREATE VEGEATION MANAGEMENT AND OBSREVATIONAL LAI TIME SERIES 
"""

fusionDict = LAICalc(aoiFile = f'{workdir}/RRNW/{site}.geojson',
					 pixelRes = S1S2BBoxSize, # size in meters of the boxes that RF will use 
					 S1FilesPath = '/Volumes/WD Elements/RRNW/VVVH',
					 S2FilesPath = '/Volumes/WD Elements/RRNW/LAI',
					 metData = met_DF[['vpd']].reset_index(),
					 firstDate = firstDate,
					 lastDate = lastDate)

"""
LOAD NEE OBSERVATIONS/MEASUREMENTS 
"""

# NEEData = pd.read_csv('/Users/vm/Desktop/dir.nosync/RRNW/RRNW_Farm_CO2Flux_Data.csv')
NEEData = pd.read_csv('/Users/vm/Desktop/dir.nosync/RRNW/ground_data/RRNW_Farm_NEE_Data.csv')
NEEData.index = pd.to_datetime(NEEData.Datetime)
# NEEData = NEEData[~NEEData.index.duplicated()]
datesList = pd.to_datetime(fusionDict['dates'])

## Expand the datetime index so that the obs_NEE covers the same time period as the met and obs_LAI TS 
## remove outliers e.g. NEE cannot exceed +/- 100 g.C.m-2
NEEData = NEEData[NEEData.Field == site]
NEEData['RECO'] = NEEData.RECO*1.03775
NEEData['GPP']  = NEEData.GPP*1.03775
NEEData['NEE'] = NEEData.RECO - NEEData.GPP
NEEData = NEEData[['NEE','GPP','RECO']].resample('D').mean()
expDFtime = pd.DataFrame(index=pd.date_range(start=datesList[0] - timedelta(days=7), 
											 end=datesList[-1] - timedelta(days=7),
											 freq='D',tz='UTC'))

D = pd.concat([expDFtime,NEEData],axis=1)

NEEDataWeek = D['NEE'].resample('7D',label='left').mean()
NEEDataWeekSD = D['NEE'].resample('7D',label='left').mean()

GPPDataWeek = D['GPP'].resample('7D',label='left').mean()
GPPDataWeekSD = D['GPP'].resample('7D',label='left').mean()

RECODataWeek = D['RECO'].resample('7D',label='left').mean()
RECODataWeekSD = D['RECO'].resample('7D',label='left').mean()

# if site == 'GreatField' : 
#     NEEData.CO2_Tower2[(NEEData.CO2_Tower2 < -50)|(NEEData.CO2_Tower2 > 50)] = np.nan 
#     NEEDataWeek = NEEData['CO2_Tower2'].resample('7D',label='left').mean() 
#     NEEDataWeekSD = NEEData['CO2_Tower2'].resample('7D',label='left').std()
# if site == 'Burrows' : 
#     NEEData.CO2_Tower4[(NEEData.CO2_Tower4 < -50)|(NEEData.CO2_Tower4 > 50)] = np.nan 
#     NEEDataWeek = NEEData['CO2_Tower4'].resample('7D',label='left').mean()
#     NEEDataWeekSD = NEEData['CO2_Tower4'].resample('7D',label='left').std() 
# if site == 'DairyField' :    
#     NEEData.CO2_Tower9[(NEEData.CO2_Tower9 < -50)|(NEEData.CO2_Tower9 > 50)] = np.nan 
#     NEEDataWeek = NEEData['CO2_Tower9'].resample('7D',label='left').mean()
#     NEEDataWeekSD = NEEData['CO2_Tower9'].resample('7D',label='left').std() 


"""
OPEN NETCDF FILE AND PROVIDE INFO ON DIMENSIONS OF VARIABLES 
"""

### Define filename
filename = f"{workdir}/RRNW/{site}_Inputs.nc" # check fileName ending 

subprocess.call(f"rm -f {filename}",shell=True) # remove .nc if it exists...

### Create a new NetCDF file
dataset = nc.Dataset(filename, "w")

### Define dimensions
metDataLen = len(weekly_DOY[:lastDate])
fwdOutsLen  = len(fusionDict['vegmg']) # CHECK, vegmg, obsLAI and obsLAIDates have dim, which is diff from metDataLen, but same as pools and fluxes 

## inputs 
dim_inputs = dataset.createDimension("dim_inputs", metDataLen)
## observations  
dim_obs = dataset.createDimension("dim_obs", metDataLen)
## pools 
dim_pools = dataset.createDimension("dim_pools", fwdOutsLen)
##fluxes 
dim_fluxes = dataset.createDimension("dim_fluxes",fwdOutsLen)
## parameters 
dim_pars = dataset.createDimension("dim_pars", posteriorLen)


"""
CREATE NETCDF VARIABLE ENTRIES 
"""


"""
MODEL INPUTS : MET , VEGETATION MANAGEMENT , OBSERVATIONAL DATA FOR LAI AND OTHER VARIABLES 
"""

## Inputs 
in_runday = dataset.createVariable("in_RunDay", "f4", ("dim_inputs",))
in_runday.long_name = "Run Day"
in_mint = dataset.createVariable("in_MinT", "f4", ("dim_inputs",))
in_mint.long_name = "Minimum Temperature"
in_mint.units = "Celcius"
in_maxt = dataset.createVariable("in_MaxT", "f4", ("dim_inputs",))
in_maxt.long_name = "Maximum Temperature"
in_maxt.units = "Celcius"
in_rad = dataset.createVariable("in_sRad", "f4", ("dim_inputs",))
in_rad.long_name = "Radiation"
in_rad.units = "MJ.m-2"
in_co2 = dataset.createVariable("in_atmCO2", "f4", ("dim_inputs",))
in_co2.long_name = "Atmospheric CO2 concentration"
in_co2.units = "Parts per million"
in_doy = dataset.createVariable("in_DOY", "f4", ("dim_inputs",))
in_doy.long_name = "Day of Year"
in_prec = dataset.createVariable("in_Prec", "f4", ("dim_inputs",))
in_prec.long_name = "Precipitation"
in_prec.units = "kg.H2O.m-2.s-1"
in_vegmg = dataset.createVariable("in_VegMgmt", "f4", ("dim_inputs",))
in_vegmg.long_name = "Vegetation Management"
in_vegmg.units = "LAI in m2.m-2 or LSU per ha"
in_baf = dataset.createVariable("in_BAF", "f4", ("dim_inputs",))
in_baf.long_name = "Burnt Area Fraction - not used"
in_mint21 = dataset.createVariable("in_MinT21", "f4", ("dim_inputs",))
in_mint21.long_name = "21-day Average Minimum Temperature"
in_mint21.units = "Kelvin"
in_photo21 = dataset.createVariable("in_Photo21", "f4", ("dim_inputs",))
in_photo21.long_name = "21-day Average Photoperiod"
in_photo21.units = "Seconds"
in_vpd21 = dataset.createVariable("in_VPD21", "f4", ("dim_inputs",))
in_vpd21.long_name = "21-day Average Vapour Pressure Deficit"
in_vpd21.units = "Pascal"
in_fmc = dataset.createVariable("in_FMC", "f4", ("dim_inputs",))
in_fmc.long_name = "Forest Management after clearing  - not used"
in_avgT = dataset.createVariable("in_avgT", "f4", ("dim_inputs",))
in_avgT.long_name = "Average T"
in_avgT.units = "C"
in_wind = dataset.createVariable("in_wind", "f4", ("dim_inputs",))
in_wind.long_name = "Average daily wind"
in_wind.units = "m.s-1"
in_vpd = dataset.createVariable("in_VPD", "f4", ("dim_inputs",))
in_vpd.long_name = "Average Vapour Pressure Deficit"
in_vpd.units = "Pascal"
in_dates = dataset.createVariable("in_dates", "str", ("dim_inputs",))
in_dates.long_name = "Date corresponding to input data entries (YYYY-MM-DD)"
in_deltat = dataset.createVariable("in_deltat", "f4", ("dim_inputs",))
in_deltat.long_name = "time step of simulation (in days)"

### Observations : LAI 
obs_LAI = dataset.createVariable("obs_LAI", "f4", ("dim_obs",))
obs_LAI.long_name = "Leaf Area Index (mean)"
obs_LAI.units = "m2.m-2"
obs_LAI_sd = dataset.createVariable("obs_LAI_sd", "f4", ("dim_obs",))
obs_LAI_sd.long_name = "Relative SD of S1 Backscatter used to create obs LAI"
obs_LAI_sd.units = "NA"
obs_LAI_dates = dataset.createVariable("obs_LAI_dates", "str", ("dim_obs",))
obs_LAI_dates.long_name = "Date corresponding to LAI obs data entries (YYYY-MM-DD)"

### Observations : NEE 
obs_NEE = dataset.createVariable("obs_NEE", "f4", ("dim_obs",))
obs_NEE.long_name = "Net Ecoystem Exchange (mean)"
obs_NEE.units = "g.C.m-2.d-1"
obs_NEE_sd = dataset.createVariable("obs_NEE_sd", "f4", ("dim_obs",))
obs_NEE_sd.long_name = "Net Ecoystem Exchange (standard deviation)"
obs_NEE_sd.units = "g.C.m-2.d-1"
obs_NEE_dates = dataset.createVariable("obs_NEE_dates", "str", ("dim_obs",))
obs_NEE_dates.long_name = "Date corresponding to NEE obs data entries (YYYY-MM-DD)"

### Observations : GPP 
obs_GPP = dataset.createVariable("obs_GPP", "f4", ("dim_obs",))
obs_GPP.long_name = "Gross Primary Production (mean)"
obs_GPP.units = "g.C.m-2.t-1"
obs_GPP_sd = dataset.createVariable("obs_GPP_sd", "f4", ("dim_obs",))
obs_GPP_sd.long_name = "Gross Primary Production (standard deviation)"
obs_GPP_sd.units = "g.C.m-2.t-1"
obs_GPP_dates = dataset.createVariable("obs_GPP_dates", "str", ("dim_obs",))
obs_GPP_dates.long_name = "Date corresponding to GPP obs data entries (YYYY-MM-DD)"

### Observations : RECO 
obs_RECO = dataset.createVariable("obs_RECO", "f4", ("dim_obs",))
obs_RECO.long_name = "Ecosystem Respiration (mean)"
obs_RECO.units = "g.C.m-2.t-1"
obs_RECO_sd = dataset.createVariable("obs_RECO_sd", "f4", ("dim_obs",))
obs_RECO_sd.long_name = "Ecosystem Respiration (standard deviation)"
obs_RECO_sd.units = "g.C.m-2.t-1"
obs_RECO_dates = dataset.createVariable("obs_RECO_dates", "str", ("dim_obs",))
obs_RECO_dates.long_name = "Date corresponding to RECO obs data entries (YYYY-MM-DD)"


"""
MODEL OUTPUTS : POOLS AND FLUXES 
"""

### FLUXES : 41 x 2 (mean and SD) = 82 vars 

flux_GPP = dataset.createVariable("flux_gpp", "f4", ("dim_fluxes",))#1 
flux_GPP_sd = dataset.createVariable("flux_gpp_sd", "f4", ("dim_fluxes",))#1 
flux_GPP.long_name = "Gross Primary Production (mean)"
flux_GPP.units = "g.C.m-2"

flux_tempRate = dataset.createVariable("flux_tempRate", "f4", ("dim_fluxes",)) #2 
flux_tempRate_sd = dataset.createVariable("flux_tempRate_sd", "f4", ("dim_fluxes",)) #2 
flux_tempRate.long_name = "Temperature Rate (mean)"
flux_tempRate.units = "NA"

flux_autoResp = dataset.createVariable("flux_autoResp", "f4", ("dim_fluxes",)) #3 
flux_autoResp_sd = dataset.createVariable("flux_autoResp_sd", "f4", ("dim_fluxes",)) #3 
flux_autoResp.long_name = "Autotrophic Respiration (mean)"
flux_autoResp.units = "NA"

flux_alloc2Leaf = dataset.createVariable("flux_alloc2Leaf", "f4", ("dim_fluxes",))#4 
flux_alloc2Leaf_sd = dataset.createVariable("flux_alloc2Leaf_sd", "f4", ("dim_fluxes",))#4 
flux_alloc2Leaf.long_name = "C allocation to leaves (mean)"
flux_alloc2Leaf.units = "g.C.m-2"

flux_alloc2Stem = dataset.createVariable("flux_alloc2Stem", "f4", ("dim_fluxes",))#5 
flux_alloc2Stem_sd = dataset.createVariable("flux_alloc2Stem_sd", "f4", ("dim_fluxes",))#5 
flux_alloc2Stem.long_name = "C allocation to stem (mean)"
flux_alloc2Stem.units = "g.C.m-2"

flux_alloc2Root = dataset.createVariable("flux_alloc2Root", "f4", ("dim_fluxes",)) #6
flux_alloc2Root_sd = dataset.createVariable("flux_alloc2Root_sd", "f4", ("dim_fluxes",)) #6
flux_alloc2Root.long_name = "C allocation to root (mean)"
flux_alloc2Root.units = "g.C.m-2"

flux_stem2leaves = dataset.createVariable("flux_stemRelease", "f4", ("dim_fluxes",))#7
flux_stem2leaves_sd = dataset.createVariable("flux_stemRelease_sd", "f4", ("dim_fluxes",))#7
flux_stem2leaves.long_name = "C transfer from stem to leaves (mean)"
flux_stem2leaves.units = "g.C.m-2"

flux_leafLitProd = dataset.createVariable("flux_leafLitProd", "f4", ("dim_fluxes",)) #9 
flux_leafLitProd_sd = dataset.createVariable("flux_leafLitProd_sd", "f4", ("dim_fluxes",)) #9 
flux_leafLitProd.long_name = "Leaf litter production (mean)"
flux_leafLitProd.units = "g.C.m-2"

flux_hetResLit = dataset.createVariable("flux_hetResLit", "f4", ("dim_fluxes",))#11
flux_hetResLit_sd = dataset.createVariable("flux_hetResLit_sd", "f4", ("dim_fluxes",))#11
flux_hetResLit.long_name = "Heterotrophic Respiration from soil litter (mean)"
flux_hetResLit.units = "g.C.m-2"

flux_hetResSOM = dataset.createVariable("flux_hetResSOM", "f4", ("dim_fluxes",)) #12
flux_hetResSOM_sd = dataset.createVariable("flux_hetResSOM_sd", "f4", ("dim_fluxes",)) #12
flux_hetResSOM.long_name = "Heterotrophic Respiration from SOM (mean)"
flux_hetResSOM.units = "g.C.m-2"

flux_lit2SOM = dataset.createVariable("flux_lit2SOM", "f4", ("dim_fluxes",))#13
flux_lit2SOM_sd = dataset.createVariable("flux_lit2SOM_sd", "f4", ("dim_fluxes",))#13
flux_lit2SOM.long_name = "C transferred from litter to SOM pool (mean)"
flux_lit2SOM.units = "g.C.m-2"

flux_leafGrowth = dataset.createVariable("flux_leafGrowth", "f4", ("dim_fluxes",))
flux_leafGrowth_sd = dataset.createVariable("flux_leafGrowth_sd", "f4", ("dim_fluxes",))
flux_leafGrowth.long_name = "leaf growth (mean)"
flux_leafGrowth.units = "g.C.m-2"

flux_GSI = dataset.createVariable("flux_GSI", "f4", ("dim_fluxes",))
flux_GSI_sd = dataset.createVariable("flux_GSI_sd", "f4", ("dim_fluxes",))
flux_GSI.long_name = "Growing Season Index (mean)"
flux_GSI.units = "unitless"

flux_animal_manure_to_soil = dataset.createVariable("flux_animal_manure_to_soil", "f4", ("dim_fluxes",))
flux_animal_manure_to_soil_sd = dataset.createVariable("flux_animal_manure_to_soil_sd", "f4", ("dim_fluxes",))
flux_animal_manure_to_soil.long_name = "animal_manure_to_soil"
flux_animal_manure_to_soil.units = "g.C.m-2"

flux_animal_respiration = dataset.createVariable("flux_animal_respiration", "f4", ("dim_fluxes",))
flux_animal_respiration_sd = dataset.createVariable("flux_animal_respiration_sd", "f4", ("dim_fluxes",))
flux_animal_respiration.long_name = "animal_respiration"
flux_animal_respiration.units = "g.C.m-2"

flux_animal_methane = dataset.createVariable("flux_animal_methane", "f4", ("dim_fluxes",))
flux_animal_methane_sd = dataset.createVariable("flux_animal_methane_sd", "f4", ("dim_fluxes",))
flux_animal_methane.long_name = "animal_methane"
flux_animal_methane.units = "g.C.m-2"

flux_harvest = dataset.createVariable("flux_harvest", "f4", ("dim_fluxes",))
flux_harvest_sd = dataset.createVariable("flux_harvest_sd", "f4", ("dim_fluxes",))
flux_harvest.long_name = "grass harvested"
flux_harvest.units = "g.C.m-2"

flux_grazed = dataset.createVariable("flux_grazed", "f4", ("dim_fluxes",))
flux_grazed_sd = dataset.createVariable("flux_grazed_sd", "f4", ("dim_fluxes",))
flux_grazed.long_name = "grass grazed"
flux_grazed.units = "g.C.m-2"

flux_evap = dataset.createVariable("flux_evap", "f4", ("dim_fluxes",))#46
flux_evap_sd = dataset.createVariable("flux_evap_sd", "f4", ("dim_fluxes",))#46
flux_evap.long_name = "Evaporation (mean)"
flux_evap.units = "NA"

flux_transpir =dataset.createVariable("flux_transpir", "f4", ("dim_fluxes",)) #47
flux_transpir_sd =dataset.createVariable("flux_transpir_sd", "f4", ("dim_fluxes",)) #47
flux_transpir.long_name = "Transpiration (mean)"
flux_transpir.units = "NA"

flux_soilEvap =dataset.createVariable("flux_soilEvap", "f4", ("dim_fluxes",)) #48
flux_soilEvap_sd =dataset.createVariable("flux_soilEvap_sd", "f4", ("dim_fluxes",)) #48
flux_soilEvap.long_name = "Soil evaporation (mean)"
flux_soilEvap.units = "NA"

flux_wetCanopEvap = dataset.createVariable("flux_wetCanopEvap", "f4", ("dim_fluxes",)) #49
flux_wetCanopEvap_sd = dataset.createVariable("flux_wetCanopEvap_sd", "f4", ("dim_fluxes",)) #49
flux_wetCanopEvap.long_name = "Wet canopy evaporation (mean)"
flux_wetCanopEvap.units = "NA"

flux_runoff = dataset.createVariable("flux_runoff", "f4", ("dim_fluxes",)) #50
flux_runoff_sd = dataset.createVariable("flux_runoff_sd", "f4", ("dim_fluxes",)) #50
flux_runoff.long_name = "Runoff (mean)"
flux_runoff.units = "NA"

flux_underflow = dataset.createVariable("flux_underflow", "f4", ("dim_fluxes",))#51
flux_underflow_sd = dataset.createVariable("flux_underflow_sd", "f4", ("dim_fluxes",))#51
flux_underflow.long_name = "Underflow (mean)"
flux_underflow.units = "NA"


### POOLS : 6 pools x 2  (mean and SD) = 12 vars 

pool_leaf = dataset.createVariable("pool_leaf", "f4", ("dim_pools",))
pool_leaf_sd = dataset.createVariable("pool_leaf_sd", "f4", ("dim_pools",))
pool_leaf.long_name = "Plant Leaves Pool"
pool_leaf.units = "g.C.m-2"

pool_stem = dataset.createVariable("pool_stem", "f4", ("dim_pools",))
pool_stem_sd = dataset.createVariable("pool_stem_sd", "f4", ("dim_pools",))
pool_stem.long_name = "Plant Stem Pool"
pool_stem.units = "g.C.m-2"

pool_root = dataset.createVariable("pool_root", "f4", ("dim_pools",))
pool_root_sd = dataset.createVariable("pool_root_sd", "f4", ("dim_pools",))
pool_root.long_name = "Plant Roots Pool"
pool_root.units = "g.C.m-2"

pool_SOM = dataset.createVariable("pool_SOM", "f4", ("dim_pools",))
pool_SOM_sd = dataset.createVariable("pool_SOM_sd", "f4", ("dim_pools",))
pool_SOM.long_name = "Soil Organic Matter Pool"
pool_SOM.units = "g.C.m-2"

pool_litter = dataset.createVariable("pool_litter", "f4", ("dim_pools",))
pool_litter_sd = dataset.createVariable("pool_litter_sd", "f4", ("dim_pools",))
pool_litter.long_name = "Soil Litter Pool"
pool_litter.units = "g.C.m-2"

pool_SW = dataset.createVariable("pool_SW", "f4", ("dim_pools",))
pool_SW_sd = dataset.createVariable("pool_SW_sd", "f4", ("dim_pools",))
pool_SW.long_name = "Soil Water Content at 0-30cm depth"
pool_SW.units = "NA"

LAI = dataset.createVariable("LAI", "f4", ("dim_pools",))
LAI_sd = dataset.createVariable("LAI_sd", "f4", ("dim_pools",))
LAI.long_name = "Leaf Area Index"
LAI.units = "m2.m-2"

GPP = dataset.createVariable("GPP", "f4", ("dim_pools",))
GPP_sd = dataset.createVariable("GPP_sd", "f4", ("dim_pools",))
GPP.long_name = "Gross Primary Production"
GPP.units = "g.C.m-2"

NEE = dataset.createVariable("NEE", "f4", ("dim_pools",))
NEE_sd = dataset.createVariable("NEE_sd", "f4", ("dim_pools",))
NEE.long_name = "Net Ecosysem Exchange"
NEE.units = "g.C.m-2"


"""
MODEL PARAMETERS : INFO AND PRIOR RANGES
"""

priors = xr.load_dataset("/Users/vm/Documents/Scripts/DALEC_Grass/DALEC/DALEC.A3.H2.M2/ManagedGrassland/example_output/RESULTS_PROCESSED/North_Wyke_states_all.nc",engine='netcdf4')    

par_P1 = dataset.createVariable('par_P1', 'f4', ('dim_pars',))
par_P1.long_name = 'Decomposition Rate'
par_P1.units = 'NA'
par_P1.def_prior_max = np.nanmax(priors.parameters[0,:,0])
par_P1.def_prior_min = np.nanmin(priors.parameters[0,:,0])
par_P1.prior_max = 0.10 
par_P1.prior_min = 1e-3

par_P2 = dataset.createVariable('par_P2', 'f4', ('dim_pars',))
par_P2.long_name = 'Fraction of GPP respired'
par_P2.units = 'P2 units'
par_P2.def_prior_max = np.nanmax(priors.parameters[0,:,1])
par_P2.def_prior_min = np.nanmin(priors.parameters[0,:,1])
par_P2.prior_max = 0.48
par_P2.prior_min = 0.43

par_P3 = dataset.createVariable('par_P3', 'f4', ('dim_pars',))
par_P3.long_name = 'GSI sensitivity for leaf growth'
par_P3.units = 'P3 units'
par_P3.def_prior_max =  np.nanmax(priors.parameters[0,:,2])
par_P3.def_prior_min = np.nanmin(priors.parameters[0,:,2])
par_P3.prior_max = 1.50
par_P3.prior_min = 0.75

par_P4 = dataset.createVariable('par_P4', 'f4', ('dim_pars',))
par_P4.long_name = 'NPP belowground allocation par_Pameter'
par_P4.units = 'P4 units'
par_P4.def_prior_max = np.nanmax(priors.parameters[0,:,3])
par_P4.def_prior_min = np.nanmin(priors.parameters[0,:,3])
par_P4.prior_max = 0.90
par_P4.prior_min = 0.10

par_P5 = dataset.createVariable('par_P5', 'f4', ('dim_pars',))
par_P5.long_name = 'GSI maximum leaf turnover'
par_P5.units = 'P5 units'
par_P5.def_prior_max = np.nanmax(priors.parameters[0,:,4])
par_P5.def_prior_min = np.nanmin(priors.parameters[0,:,4])
par_P5.prior_max = 2.0
par_P5.prior_min = 1e-3

par_P6 = dataset.createVariable('par_P6', 'f4', ('dim_pars',))
par_P6.long_name = 'TOR of roots'
par_P6.units = 'P6 units'
par_P6.def_prior_max = np.nanmax(priors.parameters[0,:,5])
par_P6.def_prior_min = np.nanmin(priors.parameters[0,:,5])
par_P6.prior_max = 1e-1
par_P6.prior_min = 1e-3

par_P7 = dataset.createVariable('par_P7', 'f4', ('dim_pars',))
par_P7.long_name = 'TOR litter'
par_P7.units = 'P7 units'
par_P7.def_prior_max = np.nanmax(priors.parameters[0,:,6])
par_P7.def_prior_min = np.nanmin(priors.parameters[0,:,6])
par_P7.prior_max = 1e-1
par_P7.prior_min = 1e-3

par_P8 = dataset.createVariable('par_P8', 'f4', ('dim_pars',))
par_P8.long_name = 'TOR SOM'
par_P8.units = 'P8 units'
par_P8.def_prior_max = np.nanmax(priors.parameters[0,:,7])
par_P8.def_prior_min = np.nanmin(priors.parameters[0,:,7])
par_P8.prior_max = 1e-4
par_P8.prior_min = 1e-7

par_P9 = dataset.createVariable('par_P9', 'f4', ('dim_pars',))
par_P9.long_name = 'Temp factor Q10'
par_P9.units = 'P9 units'
par_P9.def_prior_max = np.nanmax(priors.parameters[0,:,8])
par_P9.def_prior_min = np.nanmin(priors.parameters[0,:,8])
par_P9.prior_max = 0.20
par_P9.prior_min = 0.01

par_P10 = dataset.createVariable('par_P10', 'f4', ('dim_pars',))
par_P10.long_name = 'GSI max labile turnover'
par_P10.units = 'P10 units'
par_P10.def_prior_max = np.nanmax(priors.parameters[0,:,9])
par_P10.def_prior_min = np.nanmin(priors.parameters[0,:,9])
par_P10.prior_max = 1.0
par_P10.prior_min = 1e-3

par_P11 = dataset.createVariable('par_P11', 'f4', ('dim_pars',))
par_P11.long_name = 'Photosynthetic N use efficiency'
par_P11.units = 'P11 units'
par_P11.def_prior_max = np.nanmax(priors.parameters[0,:,10])
par_P11.def_prior_min = np.nanmin(priors.parameters[0,:,10])
par_P11.prior_max = 100 # 100 # 45 
par_P11.prior_min = 5 # 1   # 5

par_P12 = dataset.createVariable('par_P12', 'f4', ('dim_pars',))
par_P12.long_name = 'GSI min temperature threshold'
par_P12.units = 'K'
par_P12.def_prior_max = np.nanmax(priors.parameters[0,:,11])
par_P12.def_prior_min = np.nanmin(priors.parameters[0,:,11])
par_P12.prior_max = 290
par_P12.prior_min = 230

par_P13 = dataset.createVariable('par_P13', 'f4', ('dim_pars',))
par_P13.long_name = 'GSI max temperature threshold '
par_P13.units = 'K'
par_P13.def_prior_max = np.nanmax(priors.parameters[0,:,12])
par_P13.def_prior_min = np.nanmin(priors.parameters[0,:,12])
par_P13.prior_max = 300
par_P13.prior_min = 250

par_P14 = dataset.createVariable('par_P14', 'f4', ('dim_pars',))
par_P14.long_name = 'GSI min photoperiod threshold'
par_P14.units = 'seconds'
par_P14.def_prior_max = np.nanmax(priors.parameters[0,:,13])
par_P14.def_prior_min = np.nanmin(priors.parameters[0,:,13])
par_P14.prior_max = 20000
par_P14.prior_min = 3600

par_P15 = dataset.createVariable('par_P15', 'f4', ('dim_pars',))
par_P15.long_name = 'Carbon per leaf area'
par_P15.units = 'g of C per m2 of leaf area'
par_P15.def_prior_max = np.nanmax(priors.parameters[0,:,14])
par_P15.def_prior_min = np.nanmin(priors.parameters[0,:,14])
par_P15.prior_max = 100 # 60 - 100 
par_P15.prior_min = 20

par_P16 = dataset.createVariable('par_P16', 'f4', ('dim_pars',))
par_P16.long_name = 'Initial size of stem C pool'
par_P16.units = 'g.C.m-2'
par_P16.def_prior_max = np.nanmax(priors.parameters[0,:,15])
par_P16.def_prior_min = np.nanmin(priors.parameters[0,:,15])
par_P16.prior_max = 10
par_P16.prior_min = 1

par_P17 = dataset.createVariable('par_P17', 'f4', ('dim_pars',))
par_P17.long_name = 'Initial size of leaf C pool'
par_P17.units = 'g.C.m-2'
par_P17.def_prior_max = np.nanmax(priors.parameters[0,:,16])
par_P17.def_prior_min = np.nanmin(priors.parameters[0,:,16])
par_P17.prior_max = 20
par_P17.prior_min = 1

par_P18 = dataset.createVariable('par_P18', 'f4', ('dim_pars',))
par_P18.long_name = 'Initial size of root C pool'
par_P18.units ='g.C.m-2'
par_P18.def_prior_max = np.nanmax(priors.parameters[0,:,17])
par_P18.def_prior_min = np.nanmin(priors.parameters[0,:,17])
par_P18.prior_max = 300
par_P18.prior_min = 50

par_P19 = dataset.createVariable('par_P19', 'f4', ('dim_pars',))
par_P19.long_name = 'Initial size of soil litter C pool'
par_P19.units = 'g.C.m-2'
par_P19.def_prior_max = np.nanmax(priors.parameters[0,:,18])
par_P19.def_prior_min = np.nanmin(priors.parameters[0,:,18])
par_P19.prior_max = 4000
par_P19.prior_min = 1000

par_P20 = dataset.createVariable('par_P20', 'f4', ('dim_pars',))
par_P20.long_name = 'GSI max photoperiod threshold'
par_P20.units = 'seconds'
par_P20.def_prior_max = np.nanmax(priors.parameters[0,:,19])
par_P20.def_prior_min = np.nanmin(priors.parameters[0,:,19])
par_P20.prior_max = 40000
par_P20.prior_min = 10000

par_P21 = dataset.createVariable('par_P21', 'f4', ('dim_pars',))
par_P21.long_name = 'GSI min VPD threshold '
par_P21.units = 'Pascal'
par_P21.def_prior_max = np.nanmax(priors.parameters[0,:,20])
par_P21.def_prior_min = np.nanmin(priors.parameters[0,:,20])
par_P21.prior_max = 3000
par_P21.prior_min = 100

par_P22 = dataset.createVariable('par_P22', 'f4', ('dim_pars',))
par_P22.long_name = 'GSI max VPD threshold'
par_P22.units = 'Pascal'
par_P22.def_prior_max = np.nanmax(priors.parameters[0,:,21])
par_P22.def_prior_min = np.nanmin(priors.parameters[0,:,21]) 
par_P22.prior_max = 5000
par_P22.prior_min = 1000 

par_P23 = dataset.createVariable('par_P23', 'f4', ('dim_pars',))
par_P23.long_name = 'Initial size of SOM C pool' # SoilGrids based Soil C stock 0-200cm
par_P23.units = 'g.C.m-2'
par_P23.def_prior_max = soilCValue * 1.05 
par_P23.def_prior_min = soilCValue * 0.95 
par_P23.prior_max = soilCValue * 1.05 
par_P23.prior_min = soilCValue * 0.95 

par_P24 = dataset.createVariable('par_P24', 'f4', ('dim_pars',))
par_P24.long_name = 'GSI senstivity for leaf senescence '
par_P24.units = 'P24 units'
par_P24.def_prior_max = np.nanmax(priors.parameters[0,:,23])
par_P24.def_prior_min = np.nanmin(priors.parameters[0,:,23])
par_P24.prior_max = 1.00
par_P24.prior_min = 0.96

par_P25 = dataset.createVariable('par_P25', 'f4', ('dim_pars',))
par_P25.long_name = 'GSI - growing state'
par_P25.units = 'P25 units'
par_P25.def_prior_max = np.nanmax(priors.parameters[0,:,24])
par_P25.def_prior_min = np.nanmin(priors.parameters[0,:,24])
par_P25.prior_max = 3.0
par_P25.prior_min = 0.5

par_P26 = dataset.createVariable('par_P26', 'f4', ('dim_pars',))
par_P26.long_name = 'GSI - initial GSI value'
par_P26.units = 'P26 units'
par_P26.def_prior_max = np.nanmax(priors.parameters[0,:,25])
par_P26.def_prior_min = np.nanmin(priors.parameters[0,:,25])
par_P26.prior_max = 2.0
par_P26.prior_min = 1.0

par_P27 = dataset.createVariable('par_P27', 'f4', ('dim_pars',))
par_P27.long_name = 'DM minimum limit for grazing'
par_P27.units = 'g.C.m-2'
par_P27.def_prior_max = np.nanmax(priors.parameters[0,:,26])
par_P27.def_prior_min = np.nanmin(priors.parameters[0,:,26])
par_P27.prior_max = 200
par_P27.prior_min = 40

par_P28 = dataset.createVariable('par_P28', 'f4', ('dim_pars',))
par_P28.long_name = 'DM minimum limit for cutting'
par_P28.units = 'g.C.m-2'
par_P28.def_prior_max = np.nanmax(priors.parameters[0,:,27])
par_P28.def_prior_min = np.nanmin(priors.parameters[0,:,27])
par_P28.prior_max = 200
par_P28.prior_min = 50

par_P29 = dataset.createVariable('par_P29', 'f4', ('dim_pars',))
par_P29.long_name = 'leaf-vs-stem allocation factor'
par_P29.units = 'NA'
par_P29.def_prior_max = np.nanmax(priors.parameters[0,:,28])
par_P29.def_prior_min = np.nanmin(priors.parameters[0,:,28])
par_P29.prior_max = 0.95
par_P29.prior_min = 0.05

par_P30 = dataset.createVariable('par_P30', 'f4', ('dim_pars',))
par_P30.long_name = 'critical GPP for LAI growth'
par_P30.units = 'NA'
par_P30.def_prior_max = np.nanmax(priors.parameters[0,:,29])
par_P30.def_prior_min = np.nanmin(priors.parameters[0,:,29])
par_P30.prior_max = 0.50
par_P30.prior_min = 1e-3

par_P31 = dataset.createVariable('par_P31', 'f4', ('dim_pars',))
par_P31.long_name = 'DM demand of animal (% of animal weight as fraction)'
par_P31.units = 'NA'
par_P31.def_prior_max = np.nanmax(priors.parameters[0,:,30])
par_P31.def_prior_min = np.nanmin(priors.parameters[0,:,30])
par_P31.prior_max = 0.035
par_P31.prior_min = 0.010

par_P32 = dataset.createVariable('par_P32', 'f4', ('dim_pars',))
par_P32.long_name = 'Post-grazing labile loss (fraction)'
par_P32.units = 'NA'
par_P32.def_prior_max = np.nanmax(priors.parameters[0,:,31])
par_P32.def_prior_min = np.nanmin(priors.parameters[0,:,31])
par_P32.prior_max = 0.50
par_P32.prior_min = 0.01

par_P33 = dataset.createVariable('par_P33', 'f4', ('dim_pars',))
par_P33.long_name = 'Post-cut labile loss (fraction)'
par_P33.units = 'NA'
par_P33.def_prior_max = np.nanmax(priors.parameters[0,:,32])
par_P33.def_prior_min = np.nanmin(priors.parameters[0,:,32])
par_P33.prior_max = 0.75
par_P33.prior_min = 0.25

par_P34 = dataset.createVariable('par_P34', 'f4', ('dim_pars',))
par_P34.long_name = 'Minimum removed biomass to allow grazing'
par_P34.units = 'g.C.m-2.d-1'
par_P34.def_prior_max = np.nanmax(priors.parameters[0,:,33])
par_P34.def_prior_min = np.nanmin(priors.parameters[0,:,33])
par_P34.prior_max = 7.0
par_P34.prior_min = 1.0

par_P35 = dataset.createVariable('par_P35', 'f4', ('dim_pars',))
par_P35.long_name = 'Initial soil water - fraction of Field Capacity'
par_P35.units = 'P35 units'
par_P35.def_prior_max = np.nanmax(priors.parameters[0,:,34])
par_P35.def_prior_min = np.nanmin(priors.parameters[0,:,34])
par_P35.prior_max = 1.0
par_P35.prior_min = 0.5

par_P36 = dataset.createVariable('par_P36', 'f4', ('dim_pars',))
par_P36.long_name = 'Coarse root biomass needed to reach 50  '
par_P36.units = 'g of biomass per m-2'
par_P36.def_prior_max = np.nanmax(priors.parameters[0,:,35])
par_P36.def_prior_min = np.nanmin(priors.parameters[0,:,35])
par_P36.prior_max = 250
par_P36.prior_min = 10

par_P37 = dataset.createVariable('par_P37', 'f4', ('dim_pars',))
par_P37.long_name = 'Maximum rooting depth'
par_P37.units = 'm'
par_P37.def_prior_max = np.nanmax(priors.parameters[0,:,36])
par_P37.def_prior_min = np.nanmin(priors.parameters[0,:,36]) 
par_P37.prior_max = 5.0
par_P37.prior_min = 0.35


"""
ADD DATA TO VARIABLES CREATED ABOVE 
"""

vegmgDF = pd.DataFrame(index=pd.to_datetime(fusionDict['dates']),data=fusionDict['vegmg'],columns=['vegmg'])

## Sorted in the order that DGv2 wants them
in_runday[:] = np.array(np.arange(len(weekly_tmax[:lastDate])))+1
in_mint[:] = np.array(weekly_tmin[:lastDate])
in_maxt[:] = np.array(weekly_tmax[:lastDate])
in_rad[:] = np.array(weekly_rad[:lastDate])
in_co2[:] = np.array(co2_ppm[:lastDate])
in_doy[:] =  np.array(weekly_DOY[:lastDate])
in_prec[:] =  np.array(weekly_prec[:lastDate])
in_vegmg[:] = np.array(fusionDict['vegmg'])
in_baf[:] = np.zeros(len(weekly_tmax[:lastDate])) + 9999 
in_mint21[:] = np.array(weekly_21d_minT[:lastDate])
in_photo21[:] = np.array(weekly_21d_photopd[:lastDate])
in_vpd21[:] = np.array(weekly_21d_vpd[:lastDate])
in_fmc[:] = np.zeros(len(weekly_tmax[:lastDate])) + 9999 
in_avgT[:] = np.array(weekly_tavg[:lastDate])
in_wind[:] = np.array(weekly_wind[:lastDate])
in_vpd[:] = np.array(weekly_vpd[:lastDate])

## Input info time 
in_dates[:] = np.array(weekly_in_dates)
in_deltat[:] = np.zeros(len(weekly_tmax[:lastDate])) + timestep

## Observational data
obs_LAI[:] =  np.array(fusionDict['LAI_RF_mean']) 
obs_LAI_sd[:] =  np.array(fusionDict['LAI_RF_rsd'])
obs_LAI_dates[:] = np.array(fusionDict['dates'])
obs_NEE[:] = np.array(NEEDataWeek)
obs_NEE_sd[:] = np.array(NEEDataWeekSD)
obs_NEE_dates[:] = np.array(fusionDict['dates'])
obs_GPP[:] = np.array(GPPDataWeek)
obs_GPP_sd[:] = np.array(GPPDataWeekSD)
obs_GPP_dates[:] = np.array(fusionDict['dates'])
obs_RECO[:] = np.array(RECODataWeek)
obs_RECO_sd[:] = np.array(RECODataWeekSD)
obs_RECO_dates[:] = np.array(fusionDict['dates'])



"""
ADD ATTRIBUTES/METADATA AND SAVE 
"""

#### Add global attributes / metadata
dataset.area_polygon = wkt.dumps(gpd.read_file(f'{workdir}/RRNW/{site}.geojson').geometry[0])  # Example polygon
dataset.abbreviations = "GPP: Gross Primary Production. NPP: Net Primary Production. GSI: Growing Season Index. TOR: Turnover Rate. LAI: Leaf Area Index. DM: Dry Matter. VPD: Vapour Pressure Deficit" # Example duration (ISO 8601 format)
dataset.data_period = str(f"{str(weekly_DOY.index.date[0])}:{str(weekly_DOY.index.date[-1])}")  # Example duration (ISO 8601 format)
dataset.EO_LAI_RF_R2 = fusionDict['metrics'][0] 
dataset.EO_LAI_RF_RMSE = fusionDict['metrics'][1]
dataset.simulation_site_name    = site
dataset.simulation_site_lat     = round(gpd.read_file(f'{workdir}/RRNW/{site}.geojson').geometry.centroid.y[0],2)
dataset.simulation_timesteps_no = len(weekly_tmax)
dataset.simulation_start        = np.array(np.arange(len(weekly_tmax)))[0]+1
dataset.simulation_finish       = np.array(np.arange(len(weekly_tmax)))[-1]+1
dataset.simulation_spinup_days  = 52  # in weeks .. 52 = 1 year
dataset.simulation_met_no       = 16  # number of input variables incl. not used legacy ones
dataset.simulation_pools_no     = 6   # number of pools 
dataset.simulation_fluxes_no    = 51  # number of fluxes 
dataset.simulation_pars_no      = 37  # number of model parameters 

### Close the dataset
dataset.close()
print(f"NetCDF file '{filename}' created successfully!")


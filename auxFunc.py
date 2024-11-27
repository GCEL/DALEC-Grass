import os ; os.chdir('/Users/vm/Documents/Scripts/DALEC_Grass/')
import os.path
from scipy import stats
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from soilgrids import SoilGrids
from owslib.wcs import WebCoverageService
import rasterio
import xarray as xr
import geopandas as gpd 
import numpy as np
import sys 
import pandas as pd 
import datetime 
from rasterstats import zonal_stats
from shapely.geometry import shape, Point, mapping
from netCDF4 import Dataset
import geopandas as gpd
import json
from glob import glob 
import seaborn as sns 
import matplotlib.pyplot as plt 
from sentinelhub import BBoxSplitter, OsmSplitter, TileSplitter, CustomGridSplitter, UtmZoneSplitter, UtmGridSplitter, CRS
from math import sqrt
import getpass
from hda import Client, Configuration
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn import linear_model
from scipy.signal import savgol_filter
import warnings 
import cdsapi
import sys 
#warnings.filterwarnings("ignore")


def ERA5_DL(year,path2outs,AOIPolygon=[58,-6.42,50,1.71]):

	"""
	Downloads met data needed to run DG 
	year ; path2outs ; AOIPolhygon=[maxLat,minLon,minLat,maxLon]
	"""

	c = cdsapi.Client()

	c.retrieve(    
	    'reanalysis-era5-single-levels',
	    {   
	        'product_type':'reanalysis',
	        'format': 'netcdf',
	        'variable': [
	            '2m_temperature', '2m_dewpoint_temperature', 
	            'surface_solar_radiation_downwards','total_precipitation',
	            '10m_u_component_of_wind', '10m_v_component_of_wind',
	        ],
	        'year': [
	            f'{year}', 
	        ],
	        'month': [
	            '01', '02', '03',
	            '04', '05', '06',
	            '07', '08', '09',
	            '10', '11', '12',
	        ],
	        'day': [
	            '01', '02', '03',
	            '04', '05', '06',
	            '07', '08', '09',
	            '10', '11', '12',
	            '13', '14', '15',
	            '16', '17', '18',
	            '19', '20', '21',
	            '22', '23', '24',
	            '25', '26', '27',
	            '28', '29', '30',
	            '31',
	        ],
	        'time': [
	            '00:00', '01:00', '02:00',
	            '03:00', '04:00', '05:00',
	            '06:00', '07:00', '08:00',
	            '09:00', '10:00', '11:00',
	            '12:00', '13:00', '14:00',
	            '15:00', '16:00', '17:00',
	            '18:00', '19:00', '20:00',
	            '21:00', '22:00', '23:00',
	        ],
	        'area': [
	            f'{AOIPolygon[0]}',f'{AOIPolygon[1]}',
	            f'{AOIPolygon[2]}',f'{AOIPolygon[3]}',
	        ],
	    },
	    f'{path2outs}/ERA5_{year}.nc')



def wekeoLAIDL(path2json,firstDate,lastDate,path2DL):

	"""
	Queries www.wekeo.eu database for LAI data across the given area (path2json) and period (firstDate:lastDate)
	and downloads any data to the given path2DL
	* https://help.wekeo.eu/en/articles/6751608-what-is-the-hda-api-python-client-and-how-to-use-it

	Parameters
	----------
	path2json,firstDate,lastDate,path2DL : str

	Returns
	-------
	data in the path2DL
	"""

	### Configure user's credentials without a .hdarc
	conf = Configuration(user = "vmyrg", password = 'gesfo9cosqyjbUwkog')
	hda_client = Client(config = conf)

	### Location info : Lat / Lon 
	lat = round(gpd.read_file(path2json).geometry.centroid.y[0],2)
	lon = round(gpd.read_file(path2json).geometry.centroid.x[0],2)

	buffer = 0.005

	query = {
	  "datasetId": "EO:EEA:DAT:CLMS_HRVPP_VI",
	  "boundingBoxValues": [
		{
		  "name": "bbox",
		  "bbox": [
			lon - buffer,
			lat - buffer,
			lon + buffer,
			lat + buffer
		  ]
		}
	  ],
	  "dateRangeSelectValues": [
		{
		  "name": "temporal_interval",
		  "start": "2017-01-01T00:00:00.000Z",
		  "end": "2024-01-01T00:00:00.000Z"
		}
	  ],
	  "stringChoiceValues": [
		{
		  "name": "productType",
		  "value": "LAI"
		}
	  ]
	}

	# Ask the result for the query passed in parameter
	matches = hda_client.search(query)

	# List the results
	print(matches)

	Path(path2DL).mkdir(exist_ok = True)
	matches.download(download_dir=Path(path2DL).as_posix())


	query = {
	  "datasetId": "EO:EEA:DAT:CLMS_HRVPP_VI",
	  "boundingBoxValues": [
		{
		  "name": "bbox",
		  "bbox": [
			lon - buffer,
			lat - buffer,
			lon + buffer,
			lat + buffer
		  ]
		}
	  ],
	  "dateRangeSelectValues": [
		{
		  "name": "temporal_interval",
		  "start": "2017-01-01T00:00:00.000Z",
		  "end": "2024-01-01T00:00:00.000Z"
		}
	  ],
	  "stringChoiceValues": [
		{
		  "name": "productType",
		  "value": "QFLAG2"
		}
	  ]
	}

	# Ask the result for the query passed in parameter
	matches = hda_client.search(query)

	# List the results
	print(matches)

	Path(path2DL).mkdir(exist_ok = True)
	matches.download(download_dir=Path(path2DL).as_posix())


def daylength(dayOfYear, lat):
	"""Computes the length of the day (the time between sunrise and
	sunset) given the day of the year and latitude of the location.
	* https://doi.org/10.1016/0304-3800(94)00034-F

	Parameters
	----------
	dayOfYear : int
		The day of the year. 1 corresponds to 1st of January
		and 365 to 31st December (on a non-leap year).
	lat : float
		Latitude of the location in degrees. Positive values
		for north and negative for south.

	Returns
	-------
	d : float
		Daylength in hours.
	"""

	latInRad = np.deg2rad(lat)
	declinationOfEarth = 23.45*np.sin(np.deg2rad(360.0*(283.0+dayOfYear)/365.0))
	if -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
		return 24.0
	elif -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
		return 0.0
	else:
		hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))))
		return 2.0*hourAngle/15.0


def SoilGridsSOC(workdir,path2aoijson):
	"""
	Collects organic C density per soil layer from SoilGrids database 
	and returns C to total depth using the follwing equation : 
	OCD in kg/m3 : soil organic carbon density 
	HOT in m : soil horizon thickness  
	OCD [kg/m3] = ORC [%]/100 × BLD [kg/m3] × (1-CRF[%]/100) = OCS / HOT <=> OCS = HOT * OCD 
	
	Parameters 
	----------
	workdir : str 
	aoiFile : str 
	
	Returns
	-------
	soil carbon stock g.C.m-2 at 0-200cm 
	"""
	### provide area of interest (geojson)
	poly = gpd.read_file(path2aoijson)

	minx = round(poly.geometry.bounds.minx,3)[0]
	maxx = round(poly.geometry.bounds.maxx,3)[0]
	miny = round(poly.geometry.bounds.miny,3)[0]
	maxy = round(poly.geometry.bounds.maxy,3)[0] 

	resol = 10 # spatial resolution of data 
	WW = minx
	SS = miny
	EE = maxx
	NN = maxy
	
	var = 'ocd'        
	url = "http://maps.isric.org/mapserv?map=/map/{}.map".format(var)
	wcs = WebCoverageService(url, version='1.0.0')
	mean_covs = [k for k in wcs.contents.keys() if k.find("mean") != -1]
	soil_grids = SoilGrids()

	SC1 = soil_grids.get_coverage_data(service_id=var, coverage_id=mean_covs[0], 
										west=WW, south=SS, east=EE, north=NN,                                     
										width=resol, height=resol,  
										local_file=False, output=f"{workdir}/{mean_covs[0]}L1.tif",
										crs='urn:ogc:def:crs:EPSG::4326')
										# crs='urn:ogc:def:crs:EPSG::54009')
		
	
	SC2 = soil_grids.get_coverage_data(service_id=var, coverage_id=mean_covs[1], 
										west=WW, south=SS, east=EE, north=NN,                                     
										width=resol, height=resol,  
										local_file=False, output=f"{workdir}/{mean_covs[1]}L1.tif",
										crs='urn:ogc:def:crs:EPSG::4326')
										# crs='urn:ogc:def:crs:EPSG::54009')
		 

	SC3 = soil_grids.get_coverage_data(service_id=var, coverage_id=mean_covs[2], 
										west=WW, south=SS, east=EE, north=NN,                                     
										width=resol, height=resol,  
										local_file=False, output=f"{workdir}/{mean_covs[2]}L1.tif",
										crs='urn:ogc:def:crs:EPSG::4326')
										# crs='urn:ogc:def:crs:EPSG::54009')
		
	SC4 = soil_grids.get_coverage_data(service_id=var, coverage_id=mean_covs[3], 
										west=WW, south=SS, east=EE, north=NN,                                     
										width=resol, height=resol,  
										local_file=False, output=f"{workdir}/{mean_covs[3]}L1.tif",
										crs='urn:ogc:def:crs:EPSG::4326')
										# crs='urn:ogc:def:crs:EPSG::54009')
		
				
	SC5 = soil_grids.get_coverage_data(service_id=var, coverage_id=mean_covs[4], 
										west=WW, south=SS, east=EE, north=NN,                                     
										width=resol, height=resol,  
										local_file=False, output=f"{workdir}/{mean_covs[4]}L1.tif",
										crs='urn:ogc:def:crs:EPSG::4326')
										# crs='urn:ogc:def:crs:EPSG::54009')
			

	SC6 = soil_grids.get_coverage_data(service_id=var, coverage_id=mean_covs[5], 
										west=WW, south=SS, east=EE, north=NN,                                     
										width=resol, height=resol,  
										local_file=False, output=f"{workdir}/{mean_covs[5]}L1.tif",
										crs='urn:ogc:def:crs:EPSG::4326')
										# crs='urn:ogc:def:crs:EPSG::54009')
				

	return (np.nanmean(SC1.data)*0.05 + np.nanmean(SC2.data)*0.10 + np.nanmean(SC3.data)*0.15 + np.nanmean(SC4.data)*0.30 + np.nanmean(SC5.data)*0.40 + np.nanmean(SC6.data))*1000/10



def S1S2_Fusion(aoiFile,pixelRes,S1FilesPath,S2FilesPath,metData,firstDate,lastDate) : 
	"""
	Given info on area and period of interest this function 
	1. trains a machine learning (random forest) that uses S1 Radar backscatter data DOY and VPD 
	to predict S2 LAI and
	2. Uses the RF model to prodice LAI obs and vegetation reduction TS that intputs to DG 
	Ref : https://datashare.ed.ac.uk/handle/10283/4086
	
	Parameters 
	----------
	aoiFile : str 
	pixelRes : int 
	S1FilesPath : str
	S2FilesPath : str 
	metData : pandas dataframe 
	firstDate , lastDate : str 
	
	Returns
	-------
	numpy array with observational LAI data and vegetation reduction
	"""

	""" Load site polygon and split into bounding boxes """

	with open(aoiFile) as f: js = json.load(f)
	for feature in js['features']: polygon = shape(feature['geometry'])
	bbsize = int(sqrt((polygon.area * 10e5 * 10000 / pixelRes)))
	bbox_splitter = BBoxSplitter([polygon], CRS.WGS84, (bbsize,bbsize))  # bounding box will be split into bbsize * bbsize bounding boxes

	"""	Collect S1 and S2 data per sub-field box """

	if os.path.isfile(f"{aoiFile[:-8]}_S2Data_R{pixelRes}m.csv") == False : 

		### sentinel-2 LAI data directory
		folders_S2 = pd.DataFrame()
		for file in glob(f"{S2FilesPath}/*tif") : 
			folders_S2 = folders_S2._append({'path': file, 'date': (file.split("_")[3]).split('T')[0] },ignore_index=True)
		folders_S2 = folders_S2.sort_values('date')
		folders_S2.index = pd.to_datetime(folders_S2.date)
		folders_S2 = folders_S2[firstDate[:4]:lastDate[:4]]

		### collect S2 LAI data per box 
		S2_DF = pd.DataFrame()
		for y in range(len(folders_S2)) : 
			try : 
				for i in range(len((bbox_splitter.get_bbox_list()))):
					listofzones_lai = zonal_stats(bbox_splitter.bbox_list[i].geometry.wkt,folders_S2.path[y],stats=['mean','std','count'],band=1,nodata=0)
					listofzones_lai_QC = zonal_stats(bbox_splitter.bbox_list[i].geometry.wkt,folders_S2.path[y],stats=['mean','std','count'],band=2,nodata=0)
					S2_DF = S2_DF._append({'box': int(i),
										  'date': datetime.datetime.strptime(folders_S2.date[y], '%Y%m%d'),
										  'file': folders_S2.path[y],
										  'lai': listofzones_lai[0]['mean'],
										  'lai_QC': listofzones_lai_QC[0]['mean'],
										  'lai_std': listofzones_lai[0]['std']},ignore_index=True)
			except : continue 
			print(f'S2 Data : Processed file {y} / {len(folders_S2)}')

		S2_DF = S2_DF[~(S2_DF.lai_QC>0)]
		S2_DF = S2_DF[S2_DF.lai>0]
		S2_DF.to_csv(f"{aoiFile[:-8]}_S2Data_R{pixelRes}m.csv")

	else : S2_DF = pd.read_csv(f"{aoiFile[:-8]}_S2Data_R{pixelRes}m.csv")
    
	if os.path.isfile(f"{aoiFile[:-8]}_S1Data_R{pixelRes}m.csv") == False : 

		### Sentinel-1 SAR data directory 
		folders_S1 = pd.DataFrame()
		for file in glob(f"{S1FilesPath}/*tif") : 
			folders_S1 = folders_S1._append({'path': file, 'date': (file.split("_")[1][1:-4])},ignore_index=True)
		folders_S1 = folders_S1.sort_values('date')
		folders_S1.index = pd.to_datetime(folders_S1.date)
		folders_S1 = folders_S1[firstDate[:4]:lastDate[:4]]

		### collect S1 backscatter data per box 
		S1_DF = pd.DataFrame()
		for y in range(len(folders_S1)) : 
			for i in range(len((bbox_splitter.get_bbox_list()))):
				band1 = zonal_stats(bbox_splitter.bbox_list[i].wkt,folders_S1.path[y],stats=['mean','std','count'],band=1,nodata=0)
				band2 = zonal_stats(bbox_splitter.bbox_list[i].wkt,folders_S1.path[y],stats=['mean','std','count'],band=2,nodata=0)
				S1_DF = S1_DF._append({'box': int(i),
									  'date': datetime.datetime.strptime(folders_S1.date[y], '%Y%m%d'),
									  'band1': band1[0]['mean'],
									  'band1_std': band1[0]['std'],
									  'band2': band2[0]['mean'],
									  'band2_std': band2[0]['std']},ignore_index=True)
			print(f'S1 Data : Processed file {y} / {len(folders_S1)}')

		S1_DF.to_csv(f"{aoiFile[:-8]}_S1Data_R{pixelRes}m.csv")

	else : S1_DF = pd.read_csv(f"{aoiFile[:-8]}_S1Data_R{pixelRes}m.csv")


	""" Merge S1 and S2 data per sub-field box """

	### merge S1 and S2 per box 
	S1_DF = S1_DF.set_index(['date', 'box'])  # Set date and location as the index
	S2_DF = S2_DF.set_index(['date', 'box'])  # Set date and location as the index
	joined_df = pd.merge(S1_DF, S2_DF, left_index=True, right_index=True) # S1-VV/VH and S2-LAI per site box in the same DF 

	### Keep boxes/dates with S2 LAI data 
	joined_df = joined_df.reset_index()
	joined_df['DOY'] = 0
	for i in range(len(joined_df)) : joined_df['DOY'][i] = int(pd.to_datetime(joined_df['date'][i]).strftime('%j')) # add DOY column

	### Add VPD data to S1-S2 dataframe - has VV VH LAI and DOY at this line 
	# joined_df = joined_df.reset_index()
	joined_df.index = pd.to_datetime(joined_df.date)
	joined_df = joined_df[pd.to_datetime(firstDate):pd.to_datetime(lastDate)]
	s1s2dateslist = np.sort(list(set(joined_df.date)))
	joined_df['VPD'] = np.nan
	for DD in s1s2dateslist :
		joined_df['VPD'][joined_df.date==DD] = round(metData.vpd[metData.date==pd.to_datetime(DD)].iloc[0],2)


	""" Machine learning """
	
	y = joined_df.lai 
	X = joined_df[['band1','band2','DOY','VPD']]

	### Implement ML algorithms to predict LAI --  Random Forest 
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=0)
	rf = RandomForestRegressor(n_estimators=100) # 'mse'
	rf.fit(X_train, y_train)
	r2 = round(rf.score(X_test, y_test),2)
	joined_df['RF_lai'] = rf.predict(X)
	RMSE = round(sqrt(mean_squared_error(joined_df.RF_lai,joined_df.lai)),2)
	# print(f'R2: {r2} , RMSE: {RMSE}')

	
	""" Create weekly continuous obs LAI and LAI reduction TS and Export results as dictionary """

	### Add VPD and DOY cols to the S1_DF and apply the RF model for all available S1 backscatter data 
	S1_DF = S1_DF.reset_index()
	S1_DF.date = pd.to_datetime(S1_DF.date)
	S1_DF['VPD'] = np.nan
	S1_DF['DOY'] = np.nan
	metData.index = pd.to_datetime(metData.date)
	for DD in metData[firstDate:lastDate].index :
		try : 
			S1_DF.VPD[S1_DF.date==DD] = metData.vpd[metData.date==DD][0]
			S1_DF.DOY[S1_DF.date==DD] = metData.index[metData.date==DD].dayofyear[0]
		except : continue
	S1_DF = S1_DF.dropna()
	S1_DF['RF_LAI'] = rf.predict(S1_DF[['band1','band2','DOY','VPD']])    

	### Create empty continuous daily dataframe wtih date and box multiindex
	date_time_index = pd.date_range(pd.to_datetime(firstDate), pd.to_datetime(lastDate), freq='D')
	box_index = pd.RangeIndex(1, len(set(joined_df.box)))
	multi_index = pd.MultiIndex.from_tuples([(date, box) for date in date_time_index for box in box_index],names=('date', 'box'))
	empty_DF = pd.DataFrame(index=multi_index)

	## Merge empty_DF and S1_DF 
	S1_DF.index = S1_DF[['date','box']]
	S1_DF = S1_DF.set_index(['date', 'box'])
	DF = pd.merge(left=S1_DF, right=empty_DF, how='right', on=list(empty_DF.index.names))

	## Get mean and SD of RF LAI
	DF = DF.reset_index()
	DF.index = pd.to_datetime(DF.date)
	data = pd.DataFrame(index=(DF.RF_LAI.resample('D').mean()).index,columns=['RF_LAI','RF_LAI_sd'])
	data.RF_LAI = np.array(DF.RF_LAI.resample('D').mean())
	data.RF_LAI_sd = np.array(DF.RF_LAI.resample('D').std() / DF.RF_LAI.resample('D').mean())
	
	### Interpolate and smooth RF LAI TS 
	data.RF_LAI = data.RF_LAI.astype('float')
	data.RF_LAI = data.RF_LAI.interpolate(method='linear')
	# data['RF_LAI_smooth'] = dgaussian_filter1d(data.RF_LAI,6)    
	# data['RF_LAI_smooth'] = data.RF_LAI.rolling(window=15,center=True).mean()
	data['RF_LAI_smooth'] = savgol_filter(data.RF_LAI, window_length=7, polyorder=3)
	data['RF_LAI_smooth'] = data.RF_LAI_smooth.fillna(method='ffill')    
	data = data.reindex(index=pd.date_range(pd.to_datetime(str(int(firstDate) - 10000), format='%Y%m%d'),
							  pd.to_datetime(lastDate, format='%Y%m%d'),freq='D'), fill_value=np.nan)  

	### Calculate vegetation reduction
	data['vegmg'] = np.array(data.RF_LAI_smooth.diff())
	data['vegmg'][data['vegmg']>0] = 0 
	data['vegmg'] = abs(data['vegmg'])

	### Aggregate estiamtes to weekly  
	dataLAI = data.RF_LAI_smooth.resample('7D',label='left').last()
	dataVEGMG = data.vegmg.resample('7D',label='left').sum()
	dataLAIRSD = data.RF_LAI_sd.resample('7D',label='left').mean()

	### Keep original S2 LAI obs
	joined_df.index = pd.to_datetime(joined_df.date)
	obsLAIS2 = pd.DataFrame()
	obsLAIS2['lai'] = joined_df.lai.resample('D').mean()
	obsLAIS2['laiRSD'] = joined_df.lai.resample('D').std() / joined_df.lai.resample('D').mean()
	obsLAIS2 = obsLAIS2.reindex(index=pd.date_range(pd.to_datetime(str(int(firstDate) - 10000), format='%Y%m%d'),
							  pd.to_datetime(lastDate, format='%Y%m%d'),freq='D'), fill_value=np.nan)  
	dataLAIS2 = obsLAIS2['lai'].resample('7D',label='left').last()
	dataLAIS2_RSD = obsLAIS2['laiRSD'].resample('7D',label='left').last()
	
	### Create dictionary with numpy arrays as values
	data_dict = {"LAI_RF_mean": np.array(dataLAI), 
				"LAI_RF_rsd": np.array(dataLAIRSD),
				"LAI_S2_mean": np.array(dataLAIS2),
				"LAI_S2_rsd": np.array(dataLAIS2_RSD),
				"vegmg": np.array(dataVEGMG),
				"dates": np.array(dataLAI.index.strftime('%Y-%m-%d')),
				"metrics": np.array([r2, RMSE])}

	return data_dict



def figPlot(DF,site):
	"""
	Creates a figure with subplots showing obs/sim NEE LAI 
	other model outputs and inputs/drivers
	
	Parameters
	----------
	DF : pandas dataframe 
	site : str 
	
	Returns
	-------
	matplotlib figure 
	"""
	
	fig, axs = plt.subplots(6,2,figsize=(14,18))
	fig.suptitle(f'Site: {site}',fontsize=14)
	### NEE TS 
	# axs[0,0].fill_between(DF.index,DF.NEE,0)
	axs[0,0].plot(DF.NEE,color='blue')
	axs[0,0].plot(DF.NEE_map,color='green',alpha=0.5,linestyle='dashed')
	axs[0,0].scatter(DF.index,DF.obs_NEE,marker='x',color='red')
	axs[0,0].set_title('NEE')
	axs[0,0].set_ylabel('g C m-2')
	axs[0,0].legend(loc="upper right",ncol=1)
	axins0 = inset_axes(axs[0,0], width='20%', height='30%',loc='upper left',borderpad=0.5)
	axins0.scatter(DF.NEE,DF.obs_NEE,alpha=0.75,color='blue',marker='o')
	axins0.scatter(DF.NEE_map,DF.obs_NEE,alpha=0.66,color='green',marker='+')
	axins0.plot(np.arange(-15,15),alpha=0.33,linestyle='dashed',c='black')
	axins0.tick_params(axis='both', bottom=False, top=False, right=False, left=False, labelbottom=False, labelleft=False, labeltop=False, labelright=False)
	# axs[0,0].set_ylim(-5,5)
	### NEE Cumulative 
	axs[0,1].plot(DF.groupby(pd.Grouper(freq='Y'))['NEE'].cumsum(),label='annual')
	axs[0,1].plot(DF['NEE'].cumsum(),label='continuous')
	axs[0,1].legend(loc="lower left",ncol=2)
	axs[0,1].set_title('Cumulative NEE')
	axs[0,1].set_ylabel('g.C.m-2')
	### GPP 
	axs[1,0].plot(DF.GPP)
	axs[1,0].set_title('GPP')
	### Aboveground biomass  
	axs[1,1].fill_between(DF.index, DF.pool_stem,0,label='stem',color='orange')
	axs[1,1].fill_between(DF.index, DF.pool_leaf, DF.pool_stem,label='leaves',color='green',hatch='//')
	axs[1,1].fill_between(DF.index, - DF.pool_root,0,label='roots',color='brown',hatch='||')
	axs[1,1].fill_between(DF.index, - DF.pool_litter,0,label='litter',color='black',alpha=0.80)
	axs[1,1].legend(loc="upper center",ncol=4)
	axs[1,1].set_title('C Partitioning')
	axs[1,1].set_ylabel('g.C.m-2')
	axs[1,1].set_ylim(-350,300)
	axins1 = inset_axes(axs[1,1], width='25%', height='25%',loc='lower center')
	axins1.plot(DF.pool_SOM)
	axins1.set_title('Soil Organic Matter',fontsize=10)
	axins1.set_ylim(DF.pool_SOM.mean()*0.20,DF.pool_SOM.mean()*1.20)
	axins1.tick_params(axis='both', bottom=False, top=False, right=False, left=False, labelbottom=False, labelleft=True, labeltop=False, labelright=False)
	### Leaf Area Index TS 
	axs[2,0].plot(DF.obs_LAI,label='obs',color='red',alpha=0.5)
	axs[2,0].plot(DF.LAI,label='sim',color='blue')
	axs[2,0].plot(DF.LAI_map,label='MAP',color='green',alpha=0.66,linestyle='dashed')
	axs[2,0].legend(loc="upper right",ncol=1)
	axs[2,0].set_title('LAI')
	axs[2,0].set_ylabel('m2 / m2')
	axs[2,0].set_ylim(0,6)
	axins2 = inset_axes(axs[2,0], width='20%', height='30%',loc='upper left',borderpad=0.5)
	axins2.plot(np.arange(0,7),alpha=0.33,linestyle='dashed',c='black')
	axins2.scatter(DF.LAI,DF.obs_LAI,alpha=0.66,color='blue',marker='+')
	axins2.tick_params(axis='both', bottom=False, top=False, right=False, left=False, labelbottom=False, labelleft=False, labeltop=False, labelright=False)
	### Respiration flux_autoResp flux_hetResLit flux_hetResSOM
	axs[2,1].fill_between(DF.index, DF.flux_hetResLit + DF.flux_hetResSOM,0,label='Het Resp',color='blue')
	axs[2,1].fill_between(DF.index, DF.flux_autoResp, DF.flux_hetResLit + DF.flux_hetResSOM,label='Aut Resp',color='cyan')
	axs[2,1].legend()
	axs[2,1].set_title('Ecosystem Respiration')
	axs[2,1].set_ylabel('g.C.m-2')
	### Temperature
	axs[3,0].plot(DF['in_avgT'])
	axs[3,0].set_title('Temperature')
	axs[3,0].set_ylabel('C')
	### Precipitation 
	axs[3,1].plot(DF['in_Prec'])
	axs[3,1].set_title('Precipitation')
	axs[3,1].set_ylabel('mm')
	### VPD 
	axs[4,0].plot(DF['in_VPD'])
	axs[4,0].set_title('VPD')
	axs[4,0].set_ylabel('kPa')
	### SRAD 
	axs[4,1].plot( DF['in_sRad'])
	axs[4,1].set_title('Radiation (net)')
	axs[4,1].set_ylabel('MJ.m-2')
	### Evaporation 
	axs[5,0].plot( DF['flux_soilEvap'])
	axs[5,0].set_title('Soil Evaporation')
	axs[5,0].set_ylabel('')
	### Runoff 
	axs[5,1].plot( DF['flux_runoff'])
	axs[5,1].set_title('Runoff')
	axs[5,1].set_ylabel('')
	### Figure props  
	axs[0,0].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[0,1].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[1,0].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[1,1].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[2,0].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[2,1].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[3,0].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[3,1].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[4,0].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[4,1].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[5,0].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	axs[5,1].grid(color = 'gray', linestyle = '--', linewidth = 0.1)
	plt.tight_layout()
	# plt.show()

	# print(f"Accuracy : {round(100*(len(mergedDF[mergedDF.relLAIErr <= mergedDF.obs_LAI_sd])/len(mergedDF)),2)}%")
	# print(f"Relative RMSE : {round(100*sqrt(mean_squared_error(mergedDF.LAI,mergedDF.obs_LAI))/mergedDF.obs_LAI.mean(),2)}%")
	
	return 


def fwdRun(site,insFile,parsPostFile,parsSampleSize) : 
	"""
	Runs Dalec-Grass in forward mode
	
	Parameters
	----------
	site : str
	insFile : str
	parsPostFile : str
	parsSampleSize : int 
	
	Returns
	-------
	pandas dataframe 
	"""
	
	### Posteriors Info 
	prior_ranges = pd.read_csv(parsPostFile).dropna()
	prior_ranges = (prior_ranges.sort_values(by='like1',ascending=False))[:parsSampleSize]

	### Drivers and Setup 
	ins      = xr.open_dataset(insFile,engine='netcdf4')    
	met      = np.vstack([ins.in_RunDay,ins.in_MinT,ins.in_MaxT,ins.in_sRad,ins.in_atmCO2,
						  ins.in_DOY,ins.in_Prec,ins.in_VegMgmt,ins.in_BAF,ins.in_MinT21,
						  ins.in_Photo21,ins.in_VPD21,ins.in_FMC,ins.in_avgT, ins.in_wind, ins.in_VPD])
	deltat   = np.array(ins.in_deltat)
	lat      = ins.attrs['simulation_site_lat']
	nodays   = ins.attrs['simulation_timesteps_no'] 
	start    = ins.attrs['simulation_start']   
	finish   = ins.attrs['simulation_finish']   
	spinup   = ins.attrs['simulation_spinup_days']
	nomet    = ins.attrs['simulation_met_no']    
	nopools  = ins.attrs['simulation_pools_no']    
	nofluxes = ins.attrs['simulation_fluxes_no']    
	nopars   = ins.attrs['simulation_pars_no']    
	nodaysOut = nodays - spinup ### drop the 1st spinup year 

	### RUN 
	nosamples = len(prior_ranges)
	unc_DF = pd.DataFrame()

	### LAI GPP NEE
	lai_mult = np.zeros([nosamples,nodaysOut])
	gpp_mult = np.zeros([nosamples,nodaysOut])
	nee_mult = np.zeros([nosamples,nodaysOut])

	### FLUXES
	autoResp_mult    = np.zeros([nosamples,nodaysOut]) # 2
	respHetLit_mult  = np.zeros([nosamples,nodaysOut]) # 12
	respHetSOM_mult  = np.zeros([nosamples,nodaysOut]) # 13
	litter2SOM_mult  = np.zeros([nosamples,nodaysOut]) # 14

	animmanure_mult = np.zeros([nosamples,nodaysOut]) # 18
	animco2_mult = np.zeros([nosamples,nodaysOut])    # 19 
	animch4_mult = np.zeros([nosamples,nodaysOut])    # 20
	cut_mult = np.zeros([nosamples,nodaysOut])        # 21
	graze_mult = np.zeros([nosamples,nodaysOut])      # 22 

	evaporation_mult = np.zeros([nosamples,nodaysOut])   # 45 
	transpiration_mult = np.zeros([nosamples,nodaysOut]) # 46
	soilEvap_mult = np.zeros([nosamples,nodaysOut])      # 47 
	wetCanEvap_mult = np.zeros([nosamples,nodaysOut])    # 48 
	runoff_mult = np.zeros([nosamples,nodaysOut])        # 49
	underflow_mult = np.zeros([nosamples,nodaysOut])     # 50 

	### POOLS
	stemC_mult = np.zeros([nosamples,nodaysOut+1])     # 0 
	leafC_mult = np.zeros([nosamples,nodaysOut+1])     # 1
	rootC_mult = np.zeros([nosamples,nodaysOut+1])     # 2
	litterC_mult = np.zeros([nosamples,nodaysOut+1])   # 3
	soilC_mult = np.zeros([nosamples,nodaysOut+1])     # 4 
	soilWater_mult = np.zeros([nosamples,nodaysOut+1]) # 5

	### Forward runs with N parameter vectors.
	for i in range(0,nosamples):

		### Set parameter vector values
		pars = np.array(prior_ranges[prior_ranges.columns[1:(nopars+1)]].iloc[i],order="F")
		### Run DALEC-GRASS	
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

	###  Estimate mean/SD for considered variables and save to unc_DF
	for i in range(nodaysOut): 
		unc_DF = unc_DF._append({
				  ### LAI GPP NEE        
				  'LAI': stats.bayes_mvs(lai_mult[:,i],alpha=0.99)[0][0] , 
				  'LAI_sd':stats.bayes_mvs(lai_mult[:,i],alpha=0.99)[2][0],
				  'LAI_normality_metric': np.mean(lai_mult[:,i])/ np.median(lai_mult[:,i]) ,              
				  'GPP': stats.bayes_mvs(gpp_mult[:,i],alpha=0.99)[0][0] , 
				  'GPP_sd': stats.bayes_mvs(gpp_mult[:,i],alpha=0.99)[2][0],
				  'GPP_normality_metric':  np.mean(gpp_mult[:,i])/ np.median(gpp_mult[:,i]) ,             
				  'NEE': stats.bayes_mvs(nee_mult[:,i],alpha=0.99)[0][0] , 
				  'NEE_sd': stats.bayes_mvs(nee_mult[:,i],alpha=0.99)[2][0],
				  'NEE_normality_metric':  np.mean(nee_mult[:,i])/ np.median(nee_mult[:,i]) ,                
				  ### FLUXES
				  'flux_autoResp': stats.bayes_mvs(autoResp_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_autoResp_sd':  stats.bayes_mvs(autoResp_mult[:,i],alpha=0.99)[2][0],    
				  'flux_hetResLit': stats.bayes_mvs(respHetLit_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_hetResLit_sd':  stats.bayes_mvs(respHetLit_mult[:,i],alpha=0.99)[2][0],   
				  'flux_hetResSOM': stats.bayes_mvs(respHetSOM_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_hetResSOM_sd':  stats.bayes_mvs(respHetSOM_mult[:,i],alpha=0.99)[2][0],   
				  'flux_lit2SOM': stats.bayes_mvs(litter2SOM_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_lit2SOM_sd':  stats.bayes_mvs(litter2SOM_mult[:,i],alpha=0.99)[2][0],    
				  'flux_animal_manure_to_soil': stats.bayes_mvs(animmanure_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_animal_manure_to_soil_sd':  stats.bayes_mvs(animmanure_mult[:,i],alpha=0.99)[2][0],    
				  'flux_animal_respiration': stats.bayes_mvs(animco2_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_animal_respiration_sd':  stats.bayes_mvs(animco2_mult[:,i],alpha=0.99)[2][0],    
				  'flux_animal_methane': stats.bayes_mvs(animch4_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_animal_methane_sd':  stats.bayes_mvs(animch4_mult[:,i],alpha=0.99)[2][0],    
				  'flux_harvest': stats.bayes_mvs(cut_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_harvest_sd':  stats.bayes_mvs(cut_mult[:,i],alpha=0.99)[2][0],    
				  'flux_grazed': stats.bayes_mvs(graze_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_grazed_sd':  stats.bayes_mvs(graze_mult[:,i],alpha=0.99)[2][0],    
				  'flux_evap': stats.bayes_mvs(evaporation_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_evap_sd':  stats.bayes_mvs(evaporation_mult[:,i],alpha=0.99)[2][0],    
				  'flux_transpir': stats.bayes_mvs(transpiration_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_transpir_sd':  stats.bayes_mvs(transpiration_mult[:,i],alpha=0.99)[2][0],    
				  'flux_soilEvap': stats.bayes_mvs(soilEvap_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_soilEvap_sd':  stats.bayes_mvs(soilEvap_mult[:,i],alpha=0.99)[2][0],    
				  'flux_wetCanopEvap': stats.bayes_mvs(wetCanEvap_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_wetCanopEvap_sd':  stats.bayes_mvs(wetCanEvap_mult[:,i],alpha=0.99)[2][0],   
				  'flux_runoff': stats.bayes_mvs(runoff_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_runoff_sd':  stats.bayes_mvs(runoff_mult[:,i],alpha=0.99)[2][0],    
				  'flux_underflow': stats.bayes_mvs(underflow_mult[:,i],alpha=0.99)[0][0] , 
				  'flux_underflow_sd':  stats.bayes_mvs(underflow_mult[:,i],alpha=0.99)[2][0],    
				  ### POOLS
				  'pool_leaf': stats.bayes_mvs(leafC_mult[:,1+i],alpha=0.99)[0][0] , 
				  'pool_leaf_sd': stats.bayes_mvs(leafC_mult[:,1+i],alpha=0.99)[2][0],
				  'pool_stem': stats.bayes_mvs(stemC_mult[:,1+i],alpha=0.99)[0][0] , 
				  'pool_stem_sd': stats.bayes_mvs(stemC_mult[:,1+i],alpha=0.99)[2][0],
				  'pool_root': stats.bayes_mvs(rootC_mult[:,1+i],alpha=0.99)[0][0] , 
				  'pool_root_sd': stats.bayes_mvs(rootC_mult[:,1+i],alpha=0.99)[2][0],
				  'pool_litter': stats.bayes_mvs(litterC_mult[:,1+i],alpha=0.99)[0][0] , 
				  'pool_litter_sd': stats.bayes_mvs(litterC_mult[:,1+i],alpha=0.99)[2][0],  
				  'pool_SOM': stats.bayes_mvs(soilC_mult[:,1+i],alpha=0.99)[0][0] , 
				  'pool_SOM_sd': stats.bayes_mvs(soilC_mult[:,1+i],alpha=0.99)[2][0],
				  'pool_SW': stats.bayes_mvs(soilWater_mult[:,1+i],alpha=0.99)[0][0] , 
				  'pool_SW_sd': stats.bayes_mvs(soilWater_mult[:,1+i],alpha=0.99)[2][0],
		},ignore_index=True)

	## Run with MAP parameter vector 
	pars = np.array(prior_ranges[prior_ranges.columns[1:(nopars+1)]].iloc[0],order="F")
	lai,gpp,nee,fluxes,pools = DG.carbon_model_mod.carbon_model(start,finish,deltat,lat,met,pars,nopools,nofluxes)
	lai = lai[spinup:]
	gpp = gpp[spinup:]
	nee = nee[spinup:]
	unc_DF['LAI_map'] = lai
	unc_DF['GPP_map'] = gpp
	unc_DF['NEE_map'] = nee
	unc_DF['date'] = pd.to_datetime(ins.in_dates[52:])
	
	return unc_DF    
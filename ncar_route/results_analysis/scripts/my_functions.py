#!/usr/local/bin/python

#========================================================================
#========================================================================

def read_nc(infile, varname, dimension=-1, is_time=0):
	'''Read a variable from a netCDF file

	Input:
		input file path
		variable name
		dimension: if < 0, read in all dimensions of the variable; if >= 0, only read in the [dimension]th of the variable (index starts from 0). For example, if the first dimension of the variable is time, and if dimension=2, then only reads in the 3rd time step.
		is_time: if the desired variable is time (1 for time; 0 for not time). If it is time, return an array of datetime object

	Return:
		var: a numpy array of
	'''

	from netCDF4 import Dataset
	from netCDF4 import num2date

	nc = Dataset(infile, 'r')
	if is_time==0:  # if not time variable
		if dimension<0:
			var = nc.variables[varname][:]
		else:
			var = nc.variables[varname][dimension]
	if is_time==1:  # if time variable
		time = nc.variables[varname]
		if hasattr(time, 'calendar'):  # if time variable has 'calendar' attribute
			if dimension<0:
				var = num2date(time[:], time.units, time.calendar)
			else:
				var = num2date(time[dimension], time.units, time.calendar)
		else:  # if time variable does not have 'calendar' attribute
			if dimension<0:
				var = num2date(time[:], time.units)
			else:
				var = num2date(time[dimension], time.units)
	nc.close()
	return var

#========================================================================
#========================================================================

def get_nc_ts(infile, varname, timename):
	''' This function reads in a time series form a netCDF file
		(only suitable if the variable only has time dimension, or its other dimensions (e.g. lat and lon) are 1)

	Require:
		read_nc

	Input:
		infile: nc file path [string]
		varname: data variable name in the nc file [string]
		timename: time variable name in the nc file [string]

	Return:
		s: [pd.Series] object with index of time
	'''

	import pandas as pd

	data = read_nc(infile, varname, dimension=-1, is_time=0)
	time = read_nc(infile, timename, dimension=-1, is_time=1)
	data = data.squeeze()  # delete single-dimensions (e.g. lat=1, lon=1)

	s = pd.Series(data, index=time)
	return s

#==============================================================
#==============================================================

def get_nc_spatial_data(infile, varname, lonname, latname, dimension=-1):
	'''This functions reads spatial data from nc file

	Require:
		read_nc

	Input:
		infile: nc file path [string]
		varname: data variable name in the nc file [string]
		lonname, latname: lat/lon variable name in the nc file [string]
		dimension: if < 0, read in all dimensions of the variable; if >= 0, only read in the [dimension]th of the variable (index starts from 0). For example, if the first dimension of the variable is time, and if dimension=2, then only reads in the 3rd time step.

	Return:
		lons, lats: the same dimension as in the nc file
		data: if dimension=-1, the same dimension as in nc file
	'''

	data = read_nc(infile, varname, dimension=dimension, is_time=0)
	lons = read_nc(infile, lonname, dimension=-1, is_time=0)
	lats = read_nc(infile, latname, dimension=-1, is_time=0)
	return lons, lats, data

#==============================================================
#==============================================================

def find_value_index(value_array, value):
	'''Find the index of a value in an array (index starts from 0)
	Return:
		The first index of the value
	'''

	import numpy as np

	index = np.where(value_array==value)
	return index[0][0]

#==============================================================
#==============================================================

def convert_YYYYMMDD_to_datetime(year, month, day):
	''' Convert arrays of year, month, day to datetime objects 
	Input:
		year: an array of years
		month: an array of months 
		day: an array of days
		(The three arrays must be the same length)
	Return:
		A list of datetime objects
'''

	import numpy as np
	import datetime as dt

	# Check if the input arrays are the same length
	if len(year)!=len(month) or len(year)!=len(day) or len(month)!=len(day):
		print "Error: the length of input date arrays not the same length!"
		exit()

	n = len(year)
	date = []
	for i in range(n):
		date.append(dt.datetime(year=np.int(np.round(year[i])), month=np.int(np.round(month[i])), day=np.int(np.round(day[i]))))

	return date

#==============================================================
#==============================================================

def plot_date_format(ax, time_range=None, locator=None, time_format=None):
	''' This function formatting plots by plt.plot_date
	Input:
		ax: plotting axis
		time range: a tuple of two datetime objects indicating xlim. e.g., (dt.date(1991,1,1), dt.date(1992,12,31))

	'''

	import matplotlib.pyplot as plt
	import datetime as dt
	from matplotlib.dates import YearLocator, MonthLocator, DateFormatter

	# Plot time range
	if time_range!=None:
		plt.xlim(time_range[0], time_range[1])

	# Set time locator (interval)
	if locator!=None:
		if locator[0]=='year':
			ax.xaxis.set_major_locator(YearLocator(locator[1]))
		elif locator[0]=='month':
			ax.xaxis.set_major_locator(MonthLocator(interval=locator[1]))

	# Set time ticks format
	if time_format!=None:
		ax.xaxis.set_major_formatter(DateFormatter(time_format))

	return ax

#==============================================================
#==============================================================

def read_USGS_streamflow(file):
	'''This function reads USGS streamflow from the directly downloaded format (date and data are in the 3rd and 4th columns, respectively; data in cfs)

	Input: directly downloaded streamflow file path

	Return:
		date_array: a list of datetime object
		flow_array: an array of flow data [cfs]

	Note: returned data and flow might not be continuous if there is missing data!!!

	'''

	import numpy as np
	import datetime as dt

	f = open(file, 'r')
	date_array = []
	flow_array = []
	while 1:
		line = f.readline().rstrip("\n")  # read in one line
		if line=="":
			break
		if line.split()[0]=='USGS' and len(line.split())>=4:  # if data line, and data not missing
			date_string = line.split()[2]  # read in date string
			date = dt.datetime.strptime(date_string, "%Y-%m-%d")  # convert date to dt object
			date_array.append(date)

			flow = float(line.split()[3])  # read in flow and convert to float [cfs]
			flow_array.append(flow)

	flow_array = np.asarray(flow_array)
	return date_array, flow_array

#==============================================================
#==============================================================

def calc_annual_cumsum_water_year(time, data):
	'''This function calculates cumulative sum of data in each water year

	Input:
		time: corresponding datetime objects
		data: an array of data

	Return:
		time: the same as input
		data_cumsum: an array of annual (water year based) cumsum data; if there is missing data at a point, data_cumsum is np.nan after this point in this water year 
	'''

	import numpy as np

	# Check if data and time is of the same length
	if len(data)!=len(time):
		print 'Error: data and time are not of the same length!'
		exit()

	data_cumsum = np.empty(len(data))
	for i in range(len(data)):
		if i==0:  # if the first day of record
			data_cumsum[0] = data[0]
		elif time[i].month!=10 or time[i].day!=1:  # if not Oct 1st, the same water year
			if (time[i]-time[i-1]).days>1:  # if the current day is not the next day of the previous day, i.e., if there is missing data here
				print 'Warning: missing data exists!'
				data_cumsum[i] = np.nan
			else:  # if no missing data at this time step, calculate cumsum
					data_cumsum[i] = data_cumsum[i-1] + data[i]
		else:  # if the current day is Oct 1st, and i!=0, the next water year
			data_cumsum[i] = data[i]

	return time, data_cumsum

#==============================================================
#==============================================================

def convert_time_series_to_df(time, data, columns):
	'''This function converts datetime objects and data array to pandas dataframe object

	Input:
		time: a list of datetime objects, e.g. [dt.datetime(2011,1,1), dt.datetime(2011,1,3)]
		data: a 1-D or 2D array of corresponding data; if 2-D, should have the same number of rows as 'time' length
		columns: a list of column names, the same length as the number of columns of 'data', e.g. ['A', 'B', 'C']

	Return: a dataframe object
	'''

	import pandas as pd
	df = pd.DataFrame(data, index=time, columns=columns)
	return df

#==============================================================
#==============================================================

def calc_monthly_data(data):
	'''This function calculates monthly mean values

	Input: [DataFrame/Series] with index of time
	Return: a [DataFrame/Series] object, with monthly mean values (the same units as input data)
	'''

	import pandas as pd
	data_mon = data.resample("M", how='mean')
	return data_mon

#==============================================================
#==============================================================

def calc_ts_stats_by_group(data, by, stat):
	'''This function calculates statistics of time series data grouped by year, month, etc

	Input:
		df: a [pd.DataFrame/Series] object, with index of time
		by: string of group by, (select from 'year' or 'month')
		stat: statistics to be calculated, (select from 'mean')
		(e.g., if want to calculate monthly mean seasonality (12 values), by='month' and stat='mean')

	Return:
		A [dateframe/Series] object, with group as index (e.g. 1-12 for 'month')
	'''

	import pandas as pd

	if by=='year':
		if stat=='mean':
			data_result = data.groupby(lambda x:x.year).mean()
	elif by=='month':
		if stat=='mean':
			data_result = data.groupby(lambda x:x.month).mean()

	return data_result

#==============================================================
#==============================================================

def plot_format(ax, xtick_location=None, xtick_labels=None):
	'''This function formats plots by plt.plot

	Input:
		xtick_location: e.g. [1, 2, 3]
		xtick_labels: e.g. ['one', 'two', 'three']
	'''

	import matplotlib.pyplot as plt

	ax.set_xticks(xtick_location)
	ax.set_xticklabels(xtick_labels)

	return ax

#==============================================================
#==============================================================

def read_NRNI_streamflow(file, gauge_name):
	'''This function reads in NRNI csv file

	Input:
		file: [string] original NRNI csv file path (second col: date; first row: gauge_name)
		gauge_name: [string] gauge name matching the first row of the file

	Output:
		[dataframe] flow [cfs] data at this gauge with index of dates
	'''

	import pandas as pd

	df = pd.read_csv(file, skiprows=range(1,7), usecols=(1,gauge_name), index_col=0) #, parse_dates=True, date_parser=lambda x:dt.datetime.strptime(x, '%d%b%Y'))
	df.index = pd.to_datetime(df.index, format='%d%b%Y')  # make index datetime objects
	df.columns = ['flow']  # change column name to 'flow'

	return df

#==============================================================
#==============================================================

def select_time_range(data, start_datetime, end_datetime):
	''' This function selects out the part of data within a time range 

	Input:
		data: [dataframe/Series] data with index of datetime
		start_datetime: [dt.datetime] start time
		end_datetime: [dt.datetime] end time (NOTE: end_datetime itself is not included in the output; thus, end_datetime should be set a little bit later than the desired last time point)

	Return:
		Selected data (same object type as input)
	'''

	import datetime as dt

	start = data.index.searchsorted(start_datetime)
	end = data.index.searchsorted(end_datetime)

	data_selected = data.ix[start:end]

	return data_selected

#========================================================================
#========================================================================

def define_map_projection(projection='gall', llcrnrlat=-80, urcrnrlat=80, llcrnrlon=-180, urcrnrlon=180, resolution='i', land_color='grey', ocean_color='lightblue', lakes=True):
	'''Define projected map

	Return: the projection
	'''

	from mpl_toolkits.basemap import Basemap
	import matplotlib.pyplot as plt

	m = Basemap(projection=projection, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, resolution=resolution)
	m.drawlsmask(land_color=land_color, ocean_color=ocean_color, lakes=lakes)
	m.drawcoastlines(linewidth=0.75)
	m.drawstates(linewidth=0.5)
	m.drawcountries(linewidth=0.5)

	return m










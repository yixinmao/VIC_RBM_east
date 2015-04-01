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
		date.append(dt.date(year=np.int(np.round(year[i])), month=np.int(np.round(month[i])), day=np.int(np.round(day[i]))))

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











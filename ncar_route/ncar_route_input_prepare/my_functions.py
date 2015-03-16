#!/usr/local/bin/python

#==============================================================
#==============================================================

def read_nc(infile, varname, dimension=-1):
	'''Read a variable from a netCDF file

	Input:
		input file path
		variable name
		dimension: if < 0, read in all dimensions of the variable; if >= 0, only read in the [dimension]th of the variable (index starts from 0). For example, if the first dimension of the variable is time, and if dimension=2, then only reads in the 3rd time step.

	Return:
		var: a numpy array of 
	'''

	from netCDF4 import Dataset

	nc = Dataset(infile, 'r')
	if dimension<0:
		var = nc.variables[varname][:]
	else:
		var = nc.variables[varname][dimension]
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

def find_all_upstream_segs(outlet_ID, reachIndex, reachID, reachList, reachStart, reachCount):
	'''Find all upstream segments (including itself) of an outlet

	Input:
		outlet_ID: seg ID of the outlet
		reachIndex: segment index (1 to nseg)
		reachID: segment ID (corresponding to reachIndex)
		reachList: list of indice of all the upstream segments
		reachStart: start index in upstream reach listed in reachList
		reachCount: number of all the upstream segments for each segment

	Return:
		outlet_upstream_segs_ID: ID of all upstream segs (including the outlet itself)
		outlet_upstream_segs_index: indic of all upstream segments (including the outlet itself)
	'''

	import numpy as np

	# Find outlet segment index
	outlet_reachIndex = find_value_index(reachID, outlet_ID) + 1 # seg index of the outlet (starts from 1)

	# Find all upstream segment indice (including the outlet segment)
	outlet_reachStart = reachStart[outlet_reachIndex-1]  # start index (strats from 1) of the outlet seg listed in reachList
	outlet_reachCount = reachCount[outlet_reachIndex-1]  # number of upstream segments of the outlet seg listed in reachList
	outlet_upstream_segs_index = reachList[(outlet_reachStart-1):(outlet_reachStart-1+outlet_reachCount)]  # indic of all upstream segments (including the outlet itself)

	# Find corresponding segment ID of all upstream segments (including the outlet setment)
	outlet_upstream_segs_ID = reachID[outlet_upstream_segs_index-1]  # ID of all upstream segs (including the outlet itself)

	return outlet_upstream_segs_ID, outlet_upstream_segs_index

#==============================================================
#==============================================================

def find_all_hrus(segIndice, hru_id, upHruStart, upHruCount):
	'''Given some stream segments, find all the immediate contributing HRUs

	Input:
		segIndice: an array of all seg indice considered (index starts from 1)
		hru_id: an array of all hru id
		upHruStart: Start index (starts from 1) of immediate contributing HRUs for each segment listed in hru_id
		upHruCount: number of immediate upstream Hrus for each segment

	Return:
		hru_IDs_all: an array of IDs of all the immediate contributing HRUs
	'''

	import numpy as np

	hru_IDs_all = []
	for segInd in segIndice:  # for each segment
		hru_start_index = upHruStart[segInd-1]  # start index listed in hru_id, index starts from 1
		hru_count = upHruCount[segInd-1]  # number of immediate upstream HRUs of this segment
		for i in range(hru_count):  # append the contributing HRUs of this segment to the complete list
			hru_IDs_all.append(hru_id[hru_start_index-1+i])

	hru_IDs_all = np.asarray(hru_IDs_all)
	return hru_IDs_all












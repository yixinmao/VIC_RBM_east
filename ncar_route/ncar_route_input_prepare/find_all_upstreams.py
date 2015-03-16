#!/usr/local/bin/python

#================================================#
# This script finds all segments and corresponding HRUs upstream of a given stream outlet. The input to the script is a comple topography netCDF file (as input format to the NCAR routing model), and the segment ID of the stream outlet.

# In this script, 'reach' and 'segment' are synonym

# Note that segment index starts from 1
#================================================#

import numpy as np
import pickle
import my_functions

#================================================#
# Set parameters 
#================================================#
topo_path = './input/GF_ntopo_conus.nc'  # whole netCDF topography file (as input format to routing model)
outlet_ID = 17000585
output_dir = './output'  # output directory
seg_ID_attr_name = 'seg_id2'  # attribute name of seg ID in segment shapefile
hru_ID_attr_name = 'hru_id2'  # attribute name of HRU ID in HRU shapefile

#================================================#
# Load useful topology variables
#================================================#
reachIndex = my_functions.read_nc(topo_path, 'reachIndex')  # segment index (1 to nseg) !!!!!!!!!!! Note that seg index starts from 1 !!!!!!!!!!!!!
reachID = my_functions.read_nc(topo_path, 'reachID')  # segment ID (corresponding to reachIndex)
reachList = my_functions.read_nc(topo_path, 'reachList')  # list of indice of all the upstream segments
reachStart = my_functions.read_nc(topo_path, 'reachStart')  # start index in upstream reach listed in reachList
reachCount = my_functions.read_nc(topo_path, 'reachCount')  # number of all the upstream segments for each segment
hru_id = my_functions.read_nc(topo_path, 'hru_id')  # HRU ID
upHruStart = my_functions.read_nc(topo_path, 'upHruStart')  # Start index (starts from 1) of immediate upstream Hrus for each segment listed in hru_id
upHruCount = my_functions.read_nc(topo_path, 'upHruCount')  # number of immediate upstream Hrus for each segment


#================================================#
# Find all upstream segments ID
#================================================#
outlet_upstream_segs_ID, outlet_upstream_segs_index = my_functions.find_all_upstream_segs(outlet_ID, reachIndex, reachID, reachList, reachStart, reachCount) # ID of all upstream segs (including the outlet itself); # Indice of all upstream segs (including the outlet itself)

#================================================#
# Find all upstream HRU ID
#================================================#
outlet_upstream_hru_ID = my_functions.find_all_hrus(outlet_upstream_segs_index, hru_id, upHruStart, upHruCount)

#================================================#
# Convert results to arcpy format and save variables (results will be read by arcpy)
#================================================#
# Save segment ID results
qry_segID = '"%s" IN(' %seg_ID_attr_name
for i in range(len(outlet_upstream_segs_ID)):
	qry_segID = qry_segID + '%s' %outlet_upstream_segs_ID[i]
	if i<len(outlet_upstream_segs_ID)-1:  # if not the last one
		qry_segID = qry_segID + ','
qry_segID = qry_segID + ')'
f = open('%s/qry_segID_%s.pk' %(output_dir, outlet_ID), 'w')
pickle.dump(qry_segID, f)
f.close()

# Save HRU ID results
qry_hruID = '"%s" IN(' %hru_ID_attr_name
for i in range(len(outlet_upstream_hru_ID)):
	qry_hruID = qry_hruID + '%s' %outlet_upstream_hru_ID[i]
	if i<len(outlet_upstream_hru_ID)-1:  # if not the last one
		qry_hruID = qry_hruID + ','
qry_hruID = qry_hruID + ')'
f = open('%s/qry_hruID_%s.pk' %(output_dir, outlet_ID), 'w')
pickle.dump(qry_hruID, f)
f.close()








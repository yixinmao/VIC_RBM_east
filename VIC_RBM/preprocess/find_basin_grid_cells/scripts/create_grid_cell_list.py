#!/usr/local/anaconda/bin/python

import numpy as np
import my_functions
import sys

#=====================================================
# Parameter setting
#=====================================================
basin_mask_nc = sys.argv[1] # '/raid2/ymao/VIC_RBM_east/VIC_RBM/preprocess/find_basin_grid_cells/output/PNW/Yakima_Mabtom/Yakima_Mabtom.basin_mask.nc'  # basin mask nc file (typically output from RVIC)
output_grid_list = sys.argv[2] # '/raid2/ymao/VIC_RBM_east/VIC_RBM/preprocess/find_basin_grid_cells/output/PNW/Yakima_Mabtom/Yakima_Mabtom.latlon_list'
precision = sys.argv[3] # 5  # number of decimal points in the output latlon list

#=====================================================
# Parameter setting
#=====================================================
mask, lon_mesh, lat_mesh = my_functions.get_nc_spatial(basin_mask_nc, varname='mask', lonname='lon', latname='lat')

f = open(output_grid_list, 'w')
for i in range(np.shape(mask)[0]):
	for j in range(np.shape(mask)[1]):
		if mask[i,j]:
			f.write('%.*f %.*f\n' %(int(precision), lat_mesh[i,j], int(precision), lon_mesh[i,j]))
f.close()


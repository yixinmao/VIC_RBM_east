#!/usr/local/anaconda/bin/python

# This script estimates Leopold parameters for each individual site based on USGS flow velocity/width data

import argparse
import numpy as np
import my_functions

#===========================================================#
# Load in config parameters
#===========================================================#
parser = argparse.ArgumentParser()
parser.add_argument("--cfg", type=str,  help="config file for this script")
args = parser.parse_args()
cfg = my_functions.read_config(args.cfg)

# [INPUT]
# Path for original USGS field measurement data for all sites
usgs_fld_path = cfg['INPUT']['usgs_fld_path']

#===========================================================#
# Read field measurement data for all sites
#===========================================================#
# Load data (all data in ft and sec units)
df_fld = my_functions.read_USGS_fld_meas(usgs_fld_path, time_column=4, 
                                data_columns=[24,25,26,27], 
                                data_names=['discharge','width','area','velocity'],
                                code_column=2, code_name='USGS_code')
# Separate data for each site to separate smaller dataframes
# (into dictionary with site code as keys)
dict_df_fld = my_functions.separate_df_basedOnColumn(df_fld, 'USGS_code')

#===========================================================#
# Process flow velocity and calculate depth
#===========================================================#
for key in dict_df_fld.keys():
    df = dict_df_fld[key]
    np.isnan







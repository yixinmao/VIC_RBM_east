#!/usr/local/bin/python

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import my_functions

#===============================================#
# Set parameters 
#===============================================#
ncar_route_output_nc_path = '/raid2/ymao/VIC_RBM_east/ncar_route/model_run/output/from_2860/Yakama_17000585_1990_1995.nc'  # path of necCDF output file from the routing model
outlet_ID = 17000585  # outlet segment ID to be plotted and compared

vic_route_output_path = '/raid2/ymao/VIC_RBM_east/vic_route/model_run/output/Yakima/YAKIM.day'  #'/raid2/ymao/VIC_RBM_east/ncar_route/results_analysis/data/routed_flow_from_2860/vic_streamflow_daily_historical_Yakima_Mabtom_12508990.dat'  # path of VIC routing output file; [year] [month] [day] [routed flow (cfs)]
vic_route_header_line = 0  # number of header line in the file

output_plot_path = '/raid2/ymao/VIC_RBM_east/ncar_route/results_analysis/output_plots/Yakima_Mabtom_1994_my_run.png'  # output plot path

plot_start_date = dt.date(1993, 10, 1)  # start date shown on the plot
plot_end_date = dt.date(1994, 9, 30)  # end date shown on the plot
plot_flow_max = 500  # ylim on the plot, m3/s (this should be set to auto adjust!!!!!!)

time_locator = ('month', 3)  # time locator on the plot; 'year' for year; 'month' for month. e.g., ('month', 3) for plot one tick every 3 months

#===============================================#
# Load NCAR routing data, and processing
#===============================================#
# Load useful variables in nc output file from the NCAR routing model
reachID = my_functions.read_nc(ncar_route_output_nc_path, 'reachID')  # [nSeg]
NCAR_KW_flow_all = my_functions.read_nc(ncar_route_output_nc_path, 'routedRunoff') # NCAR routing, KW-PT method [time] [nSeg]; units: m3/s
NCAR_UH_flow_all = my_functions.read_nc(ncar_route_output_nc_path, 'VICroutedRunoff') # NCAR routing, DW-UH method [time] [nSeg]; units: m3/s
time_ncar = my_functions.read_nc(ncar_route_output_nc_path, 'time', is_time=1)

# Select out the flow results at the segment considered
outlet_index = my_functions.find_value_index(reachID, outlet_ID)
NCAR_KW_flow = NCAR_KW_flow_all[:,outlet_index]
NCAR_UH_flow = NCAR_UH_flow_all[:,outlet_index]

#===============================================#
# Load VIC routing data, and processing
#===============================================#
# Load data
VIC_flow_data = np.loadtxt(vic_route_output_path, skiprows=vic_route_header_line)

# Convert dates to datetime objects
date_vic = my_functions.convert_YYYYMMDD_to_datetime(VIC_flow_data[:,0], VIC_flow_data[:,1], VIC_flow_data[:,2])

# Convert flow data units
VIC_flow = VIC_flow_data[:,3]  # unit: cfs
VIC_flow = VIC_flow * np.power(12*25.4/1000, 3)  # convert unit to: m3/s

#===============================================#
# Plot time series
#===============================================#
fig = plt.figure(figsize=(12,6))
ax = plt.axes()
# plot time series
ax.plot_date(time_ncar, NCAR_KW_flow, 'b-', label='NCAR KW-PT')
ax.plot_date(time_ncar, NCAR_UH_flow, 'g-', label='NCAR DW-UH')
ax.plot_date(date_vic, VIC_flow, 'r--', label='VIC gridded route')
# add legend, label and title
plt.legend()
plt.ylabel('Flow (cms)', fontsize=16)
plt.title('Yakima Mabtom', fontsize=16)
# formatting
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date))
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date), locator=time_locator, time_format='%Y/%m')
plt.ylim(0, plot_flow_max)

fig.savefig(output_plot_path, format='png')







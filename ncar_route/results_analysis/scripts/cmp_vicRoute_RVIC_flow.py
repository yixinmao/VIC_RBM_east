#!/usr/local/bin/python

#===============================================#
# This script compares flow at an outlet routed by the original VIC routing model and RVIC 
#===============================================#

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import my_functions

#===============================================#
# Set parameters 
#===============================================#
vic_orig_route_output_path = '/raid2/ymao/VIC_RBM_east/vic_route/model_run/output/Yakima/YAKIM.day'  # path of original VIC routing output file; [year] [month] [day] [routed flow (cfs)]
vic_route_header_line = 0  # number of header line in the file

rvic_output_path = '/raid2/ymao/VIC_RBM_east/RVIC/model_run/output/Yakima/hist/yakima_v1.rvic.h0a.1996-01-01.nc'  # path of RVIC routing output file (array type of output; streamflow unit: m3/s)

output_plot_path = '/raid2/ymao/VIC_RBM_east/ncar_route/results_analysis/output_plots/cmp_vicRoute_RVIC_Yakima_Mabtom_1994.png'  # output plot path
output_diff_plot_path = '/raid2/ymao/VIC_RBM_east/ncar_route/results_analysis/output_plots/cmp_vicRoute_RVIC_diff_Yakima_Mabtom_1994.png'  # output plot path (difference)

plot_start_date = dt.date(1993, 10, 1)  # start date shown on the plot
plot_end_date = dt.date(1994, 9, 30)  # end date shown on the plot

time_locator = ('month', 3)  # time locator on the plot; 'year' for year; 'month' for month. e.g., ('month', 3) for plot one tick every 3 months

#===============================================#
# Load VIC routing data, and processing
#===============================================#
# Load data
VIC_flow_data = np.loadtxt(vic_orig_route_output_path, skiprows=vic_route_header_line)

# Convert dates to datetime objects
date_vic = my_functions.convert_YYYYMMDD_to_datetime(VIC_flow_data[:,0], VIC_flow_data[:,1], VIC_flow_data[:,2])

# Convert flow data units
VIC_flow = VIC_flow_data[:,3]  # unit: cfs
VIC_flow = VIC_flow * np.power(12*25.4/1000, 3)  # convert unit to: m3/s

#===============================================#
# Load RVIC routing data, and processing
#===============================================#
# Load useful variables in nc output file from the RVIC routing model
time_rvic = my_functions.read_nc(rvic_output_path, 'time', is_time=1)
flow_rvic = my_functions.read_nc(rvic_output_path, 'streamflow')  # m3/s
flow_rvic = flow_rvic[:,0]

#===============================================#
# Plot time series
#===============================================#
#=== plot two flow time series ===#
fig = plt.figure(figsize=(12,6))
ax = plt.axes()
# plot time series for RVIC flow (first, select out data to be plot; then, plot)
rvic_start_ind = (plot_start_date-time_rvic[0].date()).days # ind starts from 0
rvic_end_ind = (plot_end_date-time_rvic[0].date()).days
rvic_time_to_plot = time_rvic[rvic_start_ind:rvic_end_ind+1]
rvic_flow_to_plot = flow_rvic[rvic_start_ind:rvic_end_ind+1]
ax.plot_date(rvic_time_to_plot, rvic_flow_to_plot, 'b-', label='RVIC')
# plot time series for original vic flow (first, select out data to be plot; then, plot)
vic_start_ind = (plot_start_date-date_vic[0]).days # ind starts from 0
vic_end_ind = (plot_end_date-date_vic[0]).days
vic_date_to_plot = date_vic[vic_start_ind:vic_end_ind+1]
vic_flow_to_plot = VIC_flow[vic_start_ind:vic_end_ind+1]
ax.plot_date(vic_date_to_plot, vic_flow_to_plot, 'r--', label='Orig. VIC route')
# add legend, label and title
plt.legend()
plt.ylabel('Flow (cms)', fontsize=16)
plt.title('Yakima Mabtom', fontsize=16)
# formatting
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date))
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date), locator=time_locator, time_format='%Y/%m')

fig.savefig(output_plot_path, format='png')

#=== plot difference of flow ===#
fig = plt.figure(figsize=(12,6))
ax = plt.axes()
# plot time series of flow difference
time_diff_to_plot = rvic_time_to_plot
flow_diff_to_plot = rvic_flow_to_plot - vic_flow_to_plot  # RVIC - VIC
ax.plot_date(time_diff_to_plot, flow_diff_to_plot, 'k-', label='RVIC - VIC route')
# add legend, label and title
plt.legend()
plt.ylabel('Flow (cms)', fontsize=16)
plt.title('Yakima Mabtom', fontsize=16)
# formatting
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date))
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date), locator=time_locator, time_format='%Y/%m')


fig.savefig(output_diff_plot_path, format='png')






#!/usr/local/bin/python

#===============================================#
# This script compares flow at an outlet routed by the original VIC routing model and RVIC 
# Also compares both with observed streamflow
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

usgs_flow_path = '/raid2/ymao/VIC_RBM_east/ncar_route/results_analysis/data/USGS_streamflow/12508990'  # directly downloaded USGS streamflow data at the corresponding gauge

output_plot_basename = '/raid2/ymao/VIC_RBM_east/ncar_route/results_analysis/output_plots/cmp_vicRoute_RVIC_Yakima_Mabtom_1991_1995'  # output plot path basename (suffix will be added to different plots)

plot_start_date = dt.datetime(1990, 10, 1)  # start date shown on the plot (should be complete water years)
plot_end_date = dt.datetime(1995, 9, 30)  # end date shown on the plot

time_locator = ('year', 1)  # time locator on the plot; 'year' for year; 'month' for month. e.g., ('month', 3) for plot one tick every 3 months

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
# Load USGS streamflow data, and processing
#===============================================#
time_obs, flow_obs = my_functions.read_USGS_streamflow(usgs_flow_path) # flow unit: cfs
flow_obs = flow_obs * np.power(12*25.4/1000, 3)  # convert unit to: m3/s

#===============================================#
# Plot flow time series
#===============================================#
fig = plt.figure(figsize=(12,6))
ax = plt.axes()
# plot time series for RVIC flow (first, select out data to be plotted; then, plot)
rvic_start_ind = (plot_start_date.date()-time_rvic[0].date()).days # ind starts from 0
rvic_end_ind = (plot_end_date.date()-time_rvic[0].date()).days
rvic_time_to_plot = time_rvic[rvic_start_ind:rvic_end_ind+1]
rvic_flow_to_plot = flow_rvic[rvic_start_ind:rvic_end_ind+1]
ax.plot_date(rvic_time_to_plot, rvic_flow_to_plot, 'b-', label='RVIC')
# plot time series for original vic flow (first, select out data to be plotted; then, plot)
vic_start_ind = (plot_start_date.date()-date_vic[0].date()).days # ind starts from 0
vic_end_ind = (plot_end_date.date()-date_vic[0].date()).days
vic_date_to_plot = date_vic[vic_start_ind:vic_end_ind+1]
vic_flow_to_plot = VIC_flow[vic_start_ind:vic_end_ind+1]
ax.plot_date(vic_date_to_plot, vic_flow_to_plot, 'r--', label='Orig. VIC route')
# plot time series for observed flow (first, select out data to be plotted; then, plot)
# NOTE: observed flow might not be temporal continuous!!!
for i in range(len(time_obs)):
	if time_obs[i].year==plot_start_date.year and time_obs[i].month==plot_start_date.month and time_obs[i].day==plot_start_date.day:  # if the plot_start_date
		obs_start_ind = i
	elif time_obs[i].year==plot_end_date.year and time_obs[i].month==plot_end_date.month and time_obs[i].day==plot_end_date.day:  # if the plot_end_date
		obs_end_ind = i
		break
obs_time_to_plot = time_obs[obs_start_ind:obs_end_ind+1]
obs_flow_to_plot = flow_obs[obs_start_ind:obs_end_ind+1]
ax.plot_date(obs_time_to_plot, obs_flow_to_plot, 'k-', label='Observed')
# add legend, label and title
plt.legend()
plt.ylabel('Flow (cms)', fontsize=16)
plt.title('Yakima Mabtom', fontsize=16)
# formatting
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date))
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date), locator=time_locator, time_format='%Y/%m')

fig.savefig('%s.png' %output_plot_basename, format='png')

#===============================================#
# Plot difference of RVIC and orig. VIC route flow
#===============================================#
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

fig.savefig('%s_diff.png' %output_plot_basename, format='png')

#===============================================#
# Plot cumulative flow in each water year
#===============================================#
fig = plt.figure(figsize=(12,6))
ax = plt.axes()
# plot cumulative for RVIC flow
tt, rvic_flow_cumsum_to_plot = my_functions.calc_annual_cumsum_water_year(rvic_time_to_plot, rvic_flow_to_plot)
ax.plot_date(rvic_time_to_plot, rvic_flow_cumsum_to_plot, 'b-', label='RVIC')
# plot cumulative for original vic flow
tt, vic_flow_cumsum_to_plot = my_functions.calc_annual_cumsum_water_year(vic_date_to_plot, vic_flow_to_plot)
ax.plot_date(vic_date_to_plot, vic_flow_cumsum_to_plot, 'r--', label='Orig. VIC route')
# plot cumulative for observed
tt, obs_flow_cumsum_to_plot = my_functions.calc_annual_cumsum_water_year(obs_time_to_plot, obs_flow_to_plot)
ax.plot_date(obs_time_to_plot, obs_flow_cumsum_to_plot, 'k-', label='Observed')
# add legend, label and title
plt.legend()
plt.ylabel('Flow (cms)', fontsize=16)
plt.title('Cumulative flow, Yakima Mabtom', fontsize=16)
# formatting
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date))
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date), locator=time_locator, time_format='%Y/%m')

fig.savefig('%s_cumsum.png' %output_plot_basename, format='png')

#===============================================#
# Plot monthly mean flow
#===============================================#
fig = plt.figure(figsize=(12,6))
ax = plt.axes()
# calculate and plot monthly mean flow for RVIC
df_rvic_to_plot = my_functions.convert_time_series_to_df(rvic_time_to_plot, rvic_flow_to_plot, ['flow'])  # put data into dateframe
df_rvic_mon_to_plot = my_functions.calc_monthly_data(df_rvic_to_plot) # calculate monthly mean flow
ax.plot_date(df_rvic_mon_to_plot.index, df_rvic_mon_to_plot.flow, 'b-', label='RVIC')
# calculate and plot monthly mean flow for orig. VIC route
df_vic_to_plot = my_functions.convert_time_series_to_df(vic_date_to_plot, vic_flow_to_plot, ['flow']) 
df_vic_mon_to_plot = my_functions.calc_monthly_data(df_vic_to_plot)
ax.plot_date(df_vic_mon_to_plot.index, df_vic_mon_to_plot.flow, 'r--', label='Orig. VIC route')
# calculate and plot monthly mean flow for obs.
df_obs_to_plot = my_functions.convert_time_series_to_df(obs_time_to_plot, obs_flow_to_plot, ['flow'])
df_obs_mon_to_plot = my_functions.calc_monthly_data(df_obs_to_plot)
ax.plot_date(df_obs_mon_to_plot.index, df_obs_mon_to_plot.flow, 'k-', label='Observed')
# add legend, label and title
plt.legend()
plt.ylabel('Flow (cms)', fontsize=16)
plt.title('Monthly mean flow, Yakima Mabtom', fontsize=16)
# formatting
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date))
my_functions.plot_date_format(ax, time_range=(plot_start_date, plot_end_date), locator=time_locator, time_format='%Y/%m')

fig.savefig('%s_monthly.png' %output_plot_basename, format='png')





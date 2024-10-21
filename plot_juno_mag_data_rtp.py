import numpy
import datetime

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.dates as mdates

from read_juno_ephemeris_from_amda import *
from read_juno_mag_from_amda import *
from doy_to_ymd import *
from lobe_field_strength import load_lobe_field_strength

def load_juno_mag_data_rtp(time_start, time_end,
							file_ephem = "/Users/clouis/Documents/Data/ephemerides/Juno/from_AMDA/juno_jup_xyz_jso_2016_2025.txt",
							file_mag = "/Users/clouis/Documents/Data/JUNO/MAG_FGM/from_AMDA/Juno_mag/juno_fgm_bxbybz_jso_orbits_2016-2020.txt"):
	
	(date_mag,b_x,b_y,b_z) = juno_b_xyz_mag_from_amda(file_mag)
	b_x = b_x[(date_mag >= time_start) & (date_mag <= time_end)]
	b_y = b_y[(date_mag >= time_start) & (date_mag <= time_end)]
	b_z = b_z[(date_mag >= time_start) & (date_mag <= time_end)]
	b = numpy.sqrt(b_x**2 + b_y**2 + b_z**2)
	date_mag = date_mag[(date_mag >= time_start) & (date_mag <= time_end)]
	date_mag_float = datetime_to_float(date_mag)
	

	(date_ephem,x_coord,y_coord,z_coord) = juno_ephemeris_from_amda(file_ephem)
	x_coord = x_coord[(date_ephem >= time_start) & (date_ephem <= time_end)]
	y_coord = y_coord[(date_ephem >= time_start) & (date_ephem <= time_end)]
	z_coord = z_coord[(date_ephem >= time_start) & (date_ephem <= time_end)]
	date_ephem = date_ephem[(date_ephem >= time_start) & (date_ephem <= time_end)]
	date_ephem_float = datetime_to_float(date_ephem)


	x_coord_interp = numpy.interp(date_mag_float, date_ephem_float, x_coord)
	y_coord_interp = numpy.interp(date_mag_float, date_ephem_float, y_coord)
	z_coord_interp = numpy.interp(date_mag_float, date_ephem_float, z_coord)

	r = numpy.sqrt(x_coord_interp**2 + y_coord_interp**2 + z_coord_interp**2)
	th = numpy.arccos(z_coord_interp / r)
	phi = numpy.arctan2(y_coord_interp,x_coord_interp)

	b_r = b_x*numpy.sin(th)*numpy.cos(phi)+b_y*numpy.sin(th)*numpy.sin(phi)+b_z*numpy.cos(th)
	b_th = b_x * numpy.cos(th)*numpy.cos(phi) + b_y*numpy.cos(th)*numpy.sin(phi)-b_z*numpy.sin(th)
	b_phi = b_y*numpy.cos(phi) - b_x*numpy.sin(phi)

	
	date_plot = date_mag#[(date_mag >= time_start) & (date_mag <= time_end)]
	b_plot = b#[(date_mag >= time_start) & (date_mag <= time_end)]
	b_r_plot = b_r#[(date_mag >= time_start) & (date_mag <= time_end)]
	b_phi_plot = b_phi#[(date_mag >= time_start) & (date_mag <= time_end)]
	b_th_plot = b_th#[(date_mag >= time_start) & (date_mag <= time_end)]

	return(date_plot, b_plot, b_r_plot, b_th_plot, b_phi_plot)

def plot_juno_mag_data(date_plot, b_plot, b_r_plot, b_theta_plot,b_phi_plot,
						obsID = False, directory_path_out="/Users/clouis/Documents/Etudes/",
						two_panels = False, one_panel=False, legend=False,
						labelsize = 15,
						ymin = False, ymax = False, figsize=(15,10),
						add_khurana_lobes_b_value = False):
	
	if one_panel == True:
		fig, ax = plt.subplots(1,1,figsize = figsize)
		i_b_r = 0
		i_b_theta = 0
		i_b_phi = 0
		i_b_tot = 0
	elif two_panels == True:
		fig, ax = plt.subplots(2,1,sharex = True,sharey= False,figsize = figsize)
		i_b_r = 0
		i_b_theta = 0
		i_b_phi = 0
		i_b_tot = 1
	else:
		fig, ax = plt.subplots(4,1,sharex = True,sharey= False,figsize = figsize)
		i_b_r = 0
		i_b_theta = 1
		i_b_phi = 2
		i_b_tot = 3
	
	if one_panel == True:
		ax.plot(date_plot,b_r_plot, 'r', label = r"b$_r$")
		ax.plot(date_plot,b_theta_plot, 'y', label = r"b$_\theta$")
		ax.plot(date_plot,b_phi_plot, 'g', label = r"b$_\phi$")
		ax.plot(date_plot,b_plot, "b", label = r"|b|")
		if add_khurana_lobes_b_value == True:
			(date_blobe,blobe) = load_lobe_field_strength(date_plot[0], date_plot[-1])
			ax.plot(date_blobe,blobe, "k", label = r"|b| fit",linewidth=2.5)
	else:
		ax[i_b_r].plot(date_plot,b_r_plot, 'r', label = r"b$_r$")
		ax[i_b_theta].plot(date_plot,b_theta_plot, 'y', label = r"b$_\theta$")
		ax[i_b_phi].plot(date_plot,b_phi_plot, 'g', label = r"b$_\phi$")
		ax[i_b_tot].plot(date_plot,b_plot, "b", label = r"|b|")
		if add_khurana_lobes_b_value == True:
			(date_blobe,blobe) = load_lobe_field_strength(date_plot[0], date_plot[-1])
			ax[i_b_tot].plot(date_blobe,blobe, "k", label = r"|b| fit")
	
	time_start = date_plot[0]
	time_end = date_plot[-1]
	if one_panel == True:
		ax.set_ylabel(r"|B| (nT)", fontsize = labelsize+5, fontweight='bold')
		ax.set_xlabel(f'Time ({time_start.strftime("%Y-%j (%H:%M)")} to {time_end.strftime("%Y-%j (%H:%M)")})', fontsize = labelsize+5, fontweight='bold')
		if legend == True:
			ax.legend(fontsize = labelsize+15,  prop={'weight': 'bold'})

	elif two_panels == True:
		ax[i_b_r].set_ylabel(r"B (nT)", fontsize = labelsize+5)
		ax[i_b_tot].set_ylabel(r"|B| (nT)", fontsize = labelsize+5)
		ax[i_b_tot].set_xlabel(f'Time ({time_start.strftime("%Y-%j (%H:%M)")} to {time_end.strftime("%Y-%j (%H:%M)")})', fontsize = labelsize+5)
		if legend == True:
			ax[i_b_r:i_b_tot].legend(fontsize = labelsize+5)
	
	else:
		ax[i_b_r].set_ylabel(r"B$_r$ (nT)", fontsize = labelsize+5)
		ax[i_b_theta].set_ylabel(r"B$_\theta$ (nT)", fontsize = labelsize+5)
		ax[i_b_phi].set_ylabel(r"B$_\phi$ (nT)", fontsize = labelsize+5)
		ax[i_b_tot].set_ylabel(r"|B| (nT)", fontsize = labelsize+5)
		ax[i_b_tot].set_xlabel(f'Time ({time_start.strftime("%Y-%j (%H:%M)")} to {time_end.strftime("%Y-%j (%H:%M)")})', fontsize = labelsize+5)
		if legend == True:
			ax[i_b_r:i_b_tot].legend(fontsize = labelsize+5)
	
	
	plt.xlim(time_start,time_end)
	if ymin!=False and ymax!=True:
		plt.ylim(ymin,ymax)

	if one_panel == True:
		ax.tick_params(axis='x', labelsize=labelsize, width=4.5)
		ax.tick_params(axis='y', labelsize=labelsize, width=4.5)
		plt.setp(ax.spines.values(), linewidth=2.5)
		dateFmt = mdates.DateFormatter('%Y-%m-%d\n%H:%M')
		ax.xaxis.set_major_formatter(dateFmt)
		#ax.xaxis.set_major_locater()
		ax.xaxis.set_minor_locator(mdates.HourLocator())

		ax.tick_params(which='both', width=2.5)
		ax.tick_params(which='major', length=30)
		ax.tick_params(which='minor', length=15)
	else:
		ax[i_b_r].tick_params(axis='y', labelsize=labelsize, width=2.5)
		ax[i_b_theta].tick_params(axis='y', labelsize=labelsize, width=2.5)
		ax[i_b_phi].tick_params(axis='y', labelsize=labelsize, width=2.5)
		ax[i_b_tot].tick_params(axis='y', labelsize=labelsize, width=2.5)
		ax[i_b_tot].tick_params(axis='x', labelsize=labelsize, width=2.5)


		dateFmt = mdates.DateFormatter('%Y-%m-%d\n%H:%M')
		ax[i_b_r].xaxis.set_major_formatter(dateFmt)
		#ax[i_b_r].xaxis.set_major_locater()
		ax[i_b_r].xaxis.set_minor_locator(mdates.HourLocator())
		ax[i_b_theta].xaxis.set_major_formatter(dateFmt)
		#ax[i_b_theta].xaxis.set_major_locater()
		ax[i_b_theta].xaxis.set_minor_locator(mdates.HourLocator())
		ax[i_b_phi].xaxis.set_major_formatter(dateFmt)
		#ax[i_b_phi].xaxis.set_major_locater()
		ax[i_b_phi].xaxis.set_minor_locator(mdates.HourLocator())
		ax[i_b_tot].xaxis.set_major_formatter(dateFmt)
		#ax[i_b_tot].xaxis.set_major_locater()
		ax[i_b_tot].xaxis.set_minor_locator(mdates.HourLocator())
		plt.setp(ax[i_b_r:i_b_tot].spines.values(), linewidth=2.5)
	

	
	if obsID == False:
		obsID_name = ''
	else:
		obsID_name =f'_{obsID}'
	filename = directory_path_out+f'b_rtp_plot'+obsID_name+f'_{time_start.strftime("%Y%m%dT%H%M")}_{time_end.strftime("%Y%m%dT%H%M")}.pdf'
	
	plt.tight_layout()
	
	print(f"### saving plot ###")
	plt.savefig(filename, dpi = 500)
	plt.close()

	return

	

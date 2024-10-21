from scipy.io import readsav
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from read_juno_ephemeris_from_amda import *
#from plot_juno_mag_data_rtp import load_juno_mag_data_rtp
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter
import csv

def plotting_function(x,y,ax,color='k', label = "label"):
    ax = plt.plot(x,y, color = color, label = label)
    
    return


def load_lobe_field_strength(time_start, time_end, file_ephem = "/Users/clouis/Documents/Data/ephemeris/Juno/juno_jup_xyz_iau_2016_2025.txt"):
    
    (date_ephem,x_coord,y_coord,z_coord) = juno_ephemeris_from_amda(file_ephem)
    
    x_plot = x_coord[(date_ephem>=time_start)&(date_ephem <= time_end)]
    y_plot = y_coord[(date_ephem>=time_start)&(date_ephem <= time_end)]
    z_plot = z_coord[(date_ephem>=time_start)&(date_ephem <= time_end)]
    date_plot = date_ephem[(date_ephem>=time_start)&(date_ephem <= time_end)]
    R_plot = numpy.sqrt(x_plot**2+y_plot**2+z_plot**2)
    # Blobe in nT
    Blobe = 2900 * R_plot**(-1.37)
    return(date_plot,Blobe)

def save_lobe_field_strength(date_plot,Blobe, directory_path_out="./"):
    

    filename = "lobe_field_strength"+date_plot[0].strftime("%Y%m%d%H%M")+"-"+date_plot[-1].strftime("%Y%m%d%H%M")+".csv"
    with open(directory_path_out+filename, 'w') as file:
        writer = csv.writer(file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        header = ["### Date",	"B"]
        writer.writerow(header)
        for ielements in range(len(Blobe)):
            writer.writerow([date_plot[ielements].strftime("%Y-%m-%dT%H:%M"),
							Blobe[ielements]])
    
    return


def plot_lobe_field_strength(date_plot,Blobe):
    figsize = (15,10)
    labelsize = 15

    fig, ax = plt.subplots(1,1,figsize = figsize)
    ax.plot(date_plot, Blobe, color = 'black', linewidth=3.)
    ax.set_xlim(date_plot[0], date_plot[-1])
    ax.set_ylim(-12.721781063079835,13.115)
    ax.set_ylabel(r"|B| (nT)", fontsize = labelsize+5, fontweight='bold')
    ax.set_xlabel(f'Time ({date_plot[0].strftime("%Y-%j (%H:%M)")} to {date_plot[-1].strftime("%Y-%j (%H:%M)")})', fontsize = labelsize+5, fontweight='bold')
    ax.tick_params(axis='x', labelsize=labelsize, width=4.5)
    ax.tick_params(axis='y', labelsize=labelsize, width=4.5)
    plt.setp(ax.spines.values(), linewidth=2.5)
    dateFmt = mdates.DateFormatter('%Y-%m-%d\n%H:%M')
    ax.xaxis.set_major_formatter(dateFmt)
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=6))
    ax.yaxis.set_major_locator(MultipleLocator(5.000))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(which='both', width=2.5)
    ax.tick_params(which='major', length=30)
    ax.tick_params(which='minor', length=15)
    plt.tight_layout()
    plt.show()

    return()

# Need to put this function in a separate file, to avoid crossed call between this file and plot_juno_mag_data_rtp.py
#def detect_compression_from_mag_data_and_lobe_field_strength(time_start, time_end, file_mag_data="/Users/clouis/Documents/Data/JUNO/MAG_FGM/from_AMDA/Juno_mag/juno_fgm_bxbybz_iau_orbits_2016_2020.txt"):
    
#    (date_plot,Blobe) = lobe_field_strength(time_start, time_end)
    
#    (date_plot, b_plot, b_r_plot, b_th_plot, b_phi_plot) = load_juno_mag_data_rtp(time_start, time_end,
#                            file_ephem = "/Users/clouis/Documents/Data/JUNO/MAG_FGM/from_AMDA/ephemeris/juno_jup_xyz_iau_2016_2025.txt",
#                            file_mag = "/Users/clouis/Documents/Data/JUNO/MAG_FGM/from_AMDA/Juno_mag/juno_fgm_bxbybz_iau_orbits_2016-2020.txt")




#    return()
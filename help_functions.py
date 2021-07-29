# --- Do the common import 
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# --- Set some smart functions to have different colors of the labels with latex functionalities 
matplotlib.use('pgf')
matplotlib.rc('pgf', texsystem='pdflatex')  # from running latex -v
preamble = matplotlib.rcParams.setdefault('pgf.preamble', [])
preamble.append(r'\usepackage{color}')
preamble.append(r'\usepackage{lmodern}')

# Import adjust text
from adjustText import adjust_text

import numpy as np
import os
from os.path import abspath

# Package to identify the peaks
import scipy
#print(scipy.__version__)

# Import singnal packages 
from scipy.signal import chirp, find_peaks, peak_widths
from scipy.signal import argrelextrema

#Import color package for nice plots
import brewer2mpl

# --- Define the colors of the output plots
# Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
col = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


# Import round off functions to denote the sigma values in the labels
from decimal import localcontext, Decimal, ROUND_HALF_UP
with localcontext() as ctx:
    ctx.rounding = ROUND_HALF_UP


#-----------------------------------------------------------------
# --- Set plot settings 
#-----------------------------------------------------------------
#--- Set plot settings to match LaTeX<3
def init_plotting():
    plt.rcParams['figure.figsize'] = (14, 14)
    plt.rcParams['font.size'] = 45
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 20#plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
   # plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = False
    #plt.rcParams['legend.loc'] = 'center left'
    plt.rcParams['axes.linewidth'] = 1
    # --- Use latex-style on text in plots 
    plt.rcParams['text.usetex'] = True

    # --- Custumize the length of the labels
    plt.rcParams["legend.labelspacing"] = 0.2
    plt.rcParams["legend.handlelength"] = 1.0
    plt.rcParams["legend.borderaxespad"] = 0.01

    # --- Ignore warnings for generated plot
    plt.rcParams.update({'figure.max_open_warning': 0})
    
    plt.linewidth=7.0
    font_size_plot = 40 
    return font_size_plot
#-----------------------------------------------------------------

#-----------------------------------------------------------------
def elements(array):
    # --- Check if an array have elements
    return array.ndim and array.size
#-----------------------------------------------------------------


#-----------------------------------------------------------------
# --- Identify the peaks and their indices 
def find_peak_and_index(array,peak_height,window_size=0,last_peak_index=0,delta_peak=500):

    # --- Set the distance between the peaks = delta_peak approx gf.param['no_energy_points']*(2.5/100)

    
    # --- Check if the peaks should be located in a window of not 
    if window_size==0:
        
        # --- Identify peaks over the total array        
        peaks = find_peaks(array,height=peak_height,distance=delta_peak,prominence=0.004,width=30)

        # --- Write the peaks to np.arrays
        indices = np.asarray(peaks[0][:])

        #
    else:
        # --- Indetify peaks in a window
        peaks = find_peaks(array[int(last_peak_index-window_size):int(last_peak_index+window_size)],height=peak_height,distance=delta_peak,prominence=0.01,width=30)

        # --- Write the peaks to np.arrays
        indices = np.asarray(peaks[0][:])

        # --- Shift the indices back to the original array
        indices = indices + int(last_peak_index-window_size)
 
    # --- Identify the type of peaks
    if (elements(indices)==0):
        # --- No peaks detected. Set zero values 
        peak_values = 0
        peak_indices = last_peak_index
        #
    elif (elements(indices)==1):
        #----------------
        # --- One peak detected
        peak_values = array[indices[0]]
        peak_indices = indices[0]
        #----------------
    elif (elements(indices)==2):
        #----------------
        # --- Set the sizes of the return arrays
        peak_values  = np.zeros((2))
        peak_indices = [0, 0]

        # --- Two peaks detected save the highest as diople, the lowest as quadpole
        peak_values[0]  = np.amax(array[indices])
        peak_indices[0] = indices[np.argmax(array[indices])]

        # --- Delete the maximum peak
        new_indices = np.delete(indices,np.argmax(array[indices]))

        # --- Save the quadrupole peak 
        peak_values[1]  = np.amax(array[new_indices])
        peak_indices[1] = new_indices[np.argmax(array[new_indices])]

        #----------------
    elif (elements(indices)>=3):
        #----------------

        # --- Set the sizes of the return arrays
        peak_values  = np.zeros((3))
        peak_indices = [0, 0, 0]

        # --- Two peaks detected save the highest as diople, the lowest as quadpole
        peak_values[0]  = np.amax(array[indices])       
        peak_indices[0] = indices[np.argmax(array[indices])]
        print('  max index:',indices[np.argmax(array[indices])])

        # --- Delete the maximum peak
        new_indices = np.delete(indices,np.argmax(array[indices]))

        # --- Save the quadrupole peak 
        peak_values[1]  = np.amax(array[new_indices])
        peak_indices[1] = new_indices[np.argmax(array[new_indices])]

        # --- Delete the maximum peak
        new_indices_2 = np.delete(new_indices,np.argmax(array[new_indices]))

        # --- Save the quadrupole peak
        index = new_indices_2[0]
        peak_values[2]  = array[index]
        peak_indices[2] = new_indices_2[0]

        #----------------
    else :
        #----------------
        print('ERROR!')
 
        #----------------
    return peak_values, peak_indices
#-----------------------------------------------------------------

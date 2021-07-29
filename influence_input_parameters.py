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
from scipy.signal import find_peaks
#from scipy.signal import find_peaks
from scipy.signal import argrelextrema

#Import color package for nice plots
import brewer2mpl


# Import round off functions to denote the sigma values in the labels
from decimal import localcontext, Decimal, ROUND_HALF_UP
with localcontext() as ctx:
    ctx.rounding = ROUND_HALF_UP

# Current working directory is xxx/GranFilmPy/doc, we need to add xxx/GranFilmPy/src to the import path
sys.path.append('../granfilmpy/src/')
import GranFilmPy

# --- Uncomment if we want to print the parameters 
#gf = GranFilmPy.GranFilm()
#print(gf.param)

# --- Define the colors of the output plots
# Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
col = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

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

# --- Function to set the GranFilm standard parameters
def set_standard_parameters_gf():

    # --- Initialize all of the standard parameters for the type of simulation that will be performed

    #---Source 
    #  gf.param["energy_range"] = [2,4.5] # Changed from default to be consistent with the energy range
    gf.param["theta0"] = 45.
    gf.param["phi0"]   = 0. 

    #--- Geometry 
    gf.param["radius"] = 5.0
    gf.param["truncation_ratio"] = 0.12
    gf.param['media'] = ['air','al2o3p','au','al2o3p']
    #gf.param["media"] = ['air', 'al2o3p', 'ag','al2o3p'] # Materials used in the system
    gf.param["distribution"] = 'None' # 'None','Size_Distribution'


    # --- Interaction
    gf.param["arrangement"] = 'Lattice'
    gf.param["lattice_type"] = 'Square'
    gf.param["lattice_constant"] = 15 # To force the sample variable to be within physical limits
    gf.param["island_island_interaction"] = 'None' 

    #---Numerics
    gf.param['multipole_order'] = 16
    gf.param['multipole_position_ratio']=0.0
    gf.param['no_energy_points'] = 500    
    
    # --- Return the dictionary
    return gf 

def mark_peaks_in_plot(x,y,peak_height,col=col[0],peak_to_find='maximum'):
    # --- Idenitifies and add the positions and values of eventual peaks in an array y

    # --- Set an empty label list
    label_text = []
    
    # --- Calculate the wanted distance between peaks based on the number of points
    delta_peak = gf.param['no_energy_points']*(2.5/100)
    
    # --- Identify the type of peak to identify
    if peak_to_find=='maximum':        
        peak_indices = find_peaks(y,height=peak_height,distance=delta_peak)
        #peak_indices = argrelextrema(y,np.greater,order=1)
    else:    
        peak_indices = find_peaks(1/y[:],height=peak_height,distance=delta_peak,prominence=1)

    # --- Write the peaks indices to an array
    indices = np.asarray(peak_indices[0][:])
                                     
    # --- Take out values at the peak positions
    variable = x[indices]
    data     = y[indices]

    # --- Plot and mark interesting peaks 
    for peak in range(len(variable)):

        # --- Check that the peak is valid
        if round(data[peak],4)>peak_height: 
            # --- Mark the peaks in a open figure with an x 
            plt.plot(variable[peak],data[peak],'x',markersize=14,color=col)
            
            # --- Add the labels of the peaks 
            #for i,j in zip (variable,data):
                
            # --- Add the current label with color to the plot         
            label_text.append(plt.text(round(variable[peak],3),round(data[peak],3),r'$(\textcolor{black}{%s},\textcolor{black}{%s})$'%(str(round(variable[peak],3)),str(round(data[peak],3))),color=col))
                
    # --- Return the updated label 
    return label_text



def set_array_indices():
    # --- Set the definition of the integer indices
    global s, ps, p
    global par, perp 
    
    
    # Set the indices for the p, s and ps types
    s = 1
    p = 0
    ps = 2 
    # Set the indices for the parallel and perp types 
    par = 0
    perp = 1

    return


    
def elements(array):
    # --- Check if an array have elements
    return array.ndim and array.size



# --- Define the output folders
material_names='figures/'
figure_format = '.pdf'
distribution_parameter='influence_input_parameters/'


# --- Set the indices constants
set_array_indices()

# --- Set the peak limits 
peak_limit = 7.


# --- Initialize the GranFilm class 
gf = GranFilmPy.GranFilm() # Initializes a new simulation named system1.

# --- Set the standard parameters
gf = set_standard_parameters_gf()

print('media:',gf.param['media'])

# --- Save the path to the simulations 
material_names = material_names+distribution_parameter



# --- Initialize the size of the arrays
#     Macroscpic observables 
energy  = np.zeros((gf.param["no_energy_points"]))
R       = np.zeros((gf.param["no_energy_points"],3))
delta_R = np.zeros((gf.param["no_energy_points"],3))
T       = np.zeros((gf.param["no_energy_points"],3))
A       = np.zeros((gf.param["no_energy_points"],3))

#     Microscopic observables 
alpha_dipole   = np.zeros((gf.param["no_energy_points"],2))
first_suscept  = np.zeros((gf.param["no_energy_points"],2))



# --- Loop over several multipole orders M: 1 (dipole), 16 and 30
#---------------------------------------------------------------------
#  Plot the interesting quantities 
#---------------------------------------------------------------------
# --- Set the font size 
font_size = init_plotting()
#---------------------------------------------------------------------
plt.figure(dpi=150)
M = [2, 8, 16, 24]
linestyles = [':', '--', '-.','-']

#Initialize empty label tag
labels = []

# Set the peak limit
peak_limit = 0.3

# Loop over the number of multipoles 
for k in range(4):
    
    gf.param["multipole_order"] = M[k]
    
    gf()
    
    energy = gf.energy
    R[:,0],R[:,1],R[:,2] = gf.Reflectivity(mode='R')    
    plt.plot(energy, R[:,1], linestyles[k],linewidth=4.0,color=col[k],label='M={}'.format(M[k]))

    # --- Identify peaks
    #if k==0 or k==2:
    #    labels = labels + mark_peaks_in_plot(energy[:],R[:,1],peak_limit,col=col[k])


plt.legend(loc='best')
# --- Shift the labels
adjust_text(labels)
plt.xlabel(r"$\hbar\omega\ [\mathrm{eV}]$")
plt.ylabel(r'R$_s$')
plt.tight_layout()
plt.xlim(2.,4.5)
quantity_in_plot='R_s_multipoles' 
figure_name = material_names + quantity_in_plot + figure_format    
plt.savefig(figure_name,transparence=True, bbox_inches = 'tight')
#---------------------------------------------------------------------




# --- Loop over several type of particle-particle interactions: none, dipolar and quadrupolar interactions
# --- Set the font size 
font_size = init_plotting()
plt.figure(dpi=150)
linestyles = ['--', ':', '-.']
interaction = ['None', 'Dipole', 'Quadrupole']

# --- Reset the multipole order 
gf.param["multipole_order"] = 16

#Initialize empty label tag
labels = []

# Set the peak limit
peak_limit = 0.3

# --- Loop over the island_island interaction
for k in range(3):
    
    gf.param["island_island_interaction"] = interaction[k]
    
    gf()
    
    energy = gf.energy
    R[:,0],R[:,1],R[:,2] = gf.Reflectivity(mode='R')    
    plt.plot(energy, R[:,1], linestyles[k],linewidth=4.0,color=col[k],label='{}'.format(interaction[k]))

    # --- Identify peaks
    if k==0 or k==2:
        labels = labels + mark_peaks_in_plot(energy[:],R[:,1],peak_limit,col=col[k])


# --- Shift the labels
plt.legend(loc='best')
adjust_text(labels)
plt.xlabel(r"$\hbar\omega\ [\mathrm{eV}]$")
plt.ylabel(r'R$_s$')
plt.tight_layout()
plt.xlim(2.,4.5)
quantity_in_plot='R_s_island_island_interaction' 
figure_name = material_names + quantity_in_plot + figure_format    
plt.savefig(figure_name,transparence=True, bbox_inches = 'tight')
#---------------------------------------------------------------------


# --- Final system: dependence of the correction terms 
#---------------------------------------------------------------------
system2 = GranFilmPy.GranFilm() # Initialize a new system.


# # --- Configure the system
system2.param["radius"] = 5 # radius of the whole particle = 4 nm
system2.param["lattice_constant"] = 15 # distance between 2 particles, in nanometers
system2.param["truncation_ratio"] = 0.0001 # 0 for hemispheres
system2.param["media"] = ['air', 'al2o3p', 'ag','al2o3p'] # Materials used in the system
system2.param["theta0"] = 0

# --- Loop over several type of epsilon corrections: none, finite-size and quantum blue shift
# --- Set the font size 
font_size = init_plotting()
#---------------------------------------------------------------------
plt.figure(dpi=150)
plt.xlabel(r"$\hbar\omega\ [\mathrm{eV}]$")
plt.ylabel(r'R$_s$')

#--- Set empty labels 
labels = []
peak_limit = 0.15

# --- No correction, normal incidence
system2()
energy = system2.energy
Rp,Rs,Rps = system2.Reflectivity(mode='R')
plt.plot(energy,Rs, ':',linewidth=4.0,label=r'$\epsilon^B$',color=col[0])
# --- Identify peaks
labels = labels + mark_peaks_in_plot(energy[:],Rs,peak_limit,col=col[0])


# --- Finite-size correction
eps_corr_silver = {}
eps_corr_silver['Material'] = 'ag'
eps_corr_silver['Correction'] = 'Finite-Size'
eps_corr_silver['Plasma_Frequency'] = 9.17
eps_corr_silver['Damping_Frequency'] = 0.018
eps_corr_silver['Fermi_Velocity'] = 0.91
eps_corr_silver['Finite_Size_Constant'] = 0.6
system2.param.extra_parameters["ag"] = eps_corr_silver

# --- Run the system 
system2()
energy = system2.energy
Rp,Rs,Rps = system2.Reflectivity(mode='R')
plt.plot(energy,Rs, '-.',linewidth=4.0,label=r'$\epsilon^{FS}$',color=col[1])

# --- Identify peaks
labels = labels + mark_peaks_in_plot(energy[:],Rs,peak_limit,col=col[1])

# --- Surface correction
eps_corr_silver['Plasmon_Shift'] = -1.13
eps_corr_silver['Correction'] = 'Surface'
system2.param.extra_parameters["ag"] = eps_corr_silver

# --- Run the system 
system2()

Rp,Rs,Rps = system2.Reflectivity(mode='R')
plt.plot(energy,Rs, '--',linewidth=4.0,label=r'$\epsilon^{S}$',color=col[2])
# --- Identify peaks
labels = labels + mark_peaks_in_plot(energy[:],Rs,peak_limit,col=col[2])


# --- Finite-size + Surface correction
system2.param.extra_parameters['ag']['Correction'] = ['Finite-Size','Surface']
system2()

Rp,Rs,Rps = system2.Reflectivity(mode='R')
plt.plot(energy,Rs,'-',linewidth=4.0,label=r'$\epsilon^{FS,S}$',color=col[3])
# --- Identify peaks
labels = labels + mark_peaks_in_plot(energy[:],Rs,peak_limit,col=col[3])

# --- Shift the labels
plt.legend(loc='center right')
adjust_text(labels)
plt.tight_layout()
plt.xlim(2.,4.5)
quantity_in_plot='R_s_epsilon' 
figure_name = material_names + quantity_in_plot + figure_format    
plt.savefig(figure_name,transparence=True, bbox_inches = 'tight')
#---------------------------------------------------------------------


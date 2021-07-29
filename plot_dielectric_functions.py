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

# Import numpy and other necessary libraries 
import numpy as np
import os
from os.path import abspath
import scipy

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
    gf.param["energy_range"] = [2,4.5] # Changed from default to be consistent with the energy range
    gf.param["theta0"] = 45.
    gf.param["phi0"]   = 0. 

    #--- Geometry 
    gf.param["radius"] = 1
    gf.param["truncation_ratio"] = 0.001
    gf.param['media'] = ['air','mgo','ag','mgo']
    gf.param["distribution"] = 'None' # 'None','Size_Distribution'


    # --- Interaction
    gf.param["arrangement"] = 'Lattice'
    gf.param["lattice_type"] = 'Square'
    gf.param["lattice_constant"] = 31 # To force the sample variable to be within physical limits
    gf.param["island_island_interaction"] = 'None' 

    #---Numerics
    gf.param['multipole_order'] = 16
    gf.param['multipole_position_ratio']=0.0
    gf.param['no_energy_points'] = 300
    
    # --- Return without arguments 
    return gf 


# --- Define the output folders
material_names='f_'
#material_names='figures/'
figure_format = '.pdf'
distribution_parameter='_p_'
#distribution_parameter='plot_dielectric_functions/'
label_parameter = r'$\epsilon$'


# --- Initialize the GranFilm class 
gf = GranFilmPy.GranFilm() # Initializes a new simulation named system1.

# --- Set the standard parameters
gf = set_standard_parameters_gf()

# # --- Set the media for the dielectric functions
gf.param['media'] = ['air','sio2','ag','sio2']
# # --- Set the formatted streng for the label 
formatted_media = ['\mathrm{Air}','\mathrm{SiO}_2','\mathrm{Ag}','\mathrm{SiO}_2']


# # --- Set the media for the dielectric functions
# gf.param['media'] = ['air','mgo','ag','mgo']
# # --- Set the formatted streng for the label 
# formatted_media = ['\mathrm{Air}','\mathrm{MgO}','\mathrm{Ag}','\mathrm{MgO}']


# --- Set the media for the dielectric functions
#gf.param['media'] = ['air','al2o3','ag','al2o3']
# --- Set the formatted streng for the label 
#formatted_media = ['\mathrm{Air}','\mathrm{Al}_2\mathrm{O}_3','\mathrm{Ag}','\mathrm{Al}_2\mathrm{O}_3']



# --- Save the path to the simulations 
material_names = material_names+distribution_parameter

# --- The number of different epsilons to plot
region_number = np.size(gf.param['media'])

# --- Epsilon corrections arrays
epsilon_par          = np.zeros((region_number,gf.param["no_energy_points"]),dtype=complex)
dielectric_energy    = np.zeros((region_number,gf.param["no_energy_points"]))


#---------------------------------------------------------------------
# --- Plot the dielectric function for different materials

# --- Run GranFilm to obtain the dielectric functions 
gf()

# --- Save the interesting quantities
dielectric_energy[0,:]  = gf.energy

# --- Take out the epsilon values from the different regions inside the particle 
for i in range(0,(region_number-1)):

    print(i)
    # --- Update the epsilon inside the particle 
    epsilon_par[i,:] = gf.Epsilon(gf.param['media'][i])

    print('media:',gf.param['media'][i])


#---------------------------------------------------------------------
#  Plot the interesting quantities 
#---------------------------------------------------------------------
# --- Set the font size 
font_size = init_plotting()
#---------------------------------------------------------------------

# --- The number of points to include in the plot
n=5


# --- Take out the epsilon values from the different regions inside the particle 
for i in range(0,(region_number-1)): 
    #---------------------------------------------------------------------
    # Plot dielectric function for the surface sorrection     
    plt.figure(dpi=120)
    plt.axhline(y=0., linestyle='--',linewidth=1.0,color='gray')
    
    # --- Plot the dielectric function for the current media 
    plt.scatter(dielectric_energy[0,::n],epsilon_par[i,::n].real,marker='.',linewidth=4.0,label=r'$\mathrm{Re}$',color=col[0])
    plt.scatter(dielectric_energy[0,::n],epsilon_par[i,::n].imag,marker=',',linewidth=4.0,label=r'$\mathrm{Im}$',color=col[1])
    plt.xlabel(r"$\hbar\omega\ [\mathrm{eV}]$")
    plt.ylabel(r'$\epsilon_{{{}}}$'.format(formatted_media[i]))
    plt.legend(loc='best')
    plt.xlim(2,4.5)
    quantity_in_plot_comma='epsilon_'+gf.param['media'][i]
    quantity_in_plot=quantity_in_plot_comma.replace('.','_')                        
    figure_name = material_names + quantity_in_plot + figure_format
    plt.savefig(figure_name, dpi=200, transparence=True, bbox_inches = 'tight')
    #---------------------------------------------------------------------

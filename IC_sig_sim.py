# -*- coding: utf-8 -*-
"""
Michał Włodarczyk

Script automating simulations of the signals for the Ionisation Chamber with SRIM.
Includes methods for error estimation and tools to predict if particular beam
components can be distinguished with the IC (based on the empirically found 
resolution; note that due to errors in the setup of the detector during the 2025
experimental campaign, the resolution may be slightly overestimated, 
so the results should be treated as a conservative estimate).

Pysrim package is a bit outdated and in order for the code to run correctly, one may need to find
.../site-packages/srim/core/elementdb.py and change line "yaml.load(open(dbpath, "r"))" to 
"yaml.safe_load(open(dbpath, "r"))" (the former is deprecated)

"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from srim import TRIM, Ion, Layer, Target #requires the srim package
from nuclyr import mass
from periodictable import elements
from scipy.integrate import quad
import os

import configparser #convenient package for configfiles
import collections
import collections.abc
collections.Hashable = collections.abc.Hashable #THIS IS NECESSARY FOR BACKWARDS COMPATIBILITY

atomic_mass_unit =  931.49410372 #atomic mass unit in MeV


def get_settings(settings_filename):
    """
    Imports the settings from settings_filename and stores them inside an object
    from the configparser class

    Parameters
    ----------
    settings_filename : string

    Returns
    -------
    settings : configparser

    """
    settings = configparser.ConfigParser(inline_comment_prefixes="#")
    settings.read(settings_filename)
    return settings

def gas_density_calc(molar_mass, pressure, temp = 293, scaling_factor = 1) :
    """
    Get the density of a gas given its molar mass and pressure. Uses the ideal gas law
    and assumes room temperature unless specified otherwise. (using the Van Der Waals equation
    was shown to have no major effects; for CF4 gas at least)

    Parameters
    ----------
    molar_mass : float
    pressure : float
        In mbar.
    temp : float
        In Kelvin. Defaults to room temperature.
    Returns
    -------
    density : float
        In g/cm^3 (preffered units of SRIM).

    """
    gas_const = 8.31
    return scaling_factor * molar_mass * pressure / (temp * gas_const) / 10**4

def define_target(settings, simulation_uncertainty = 1):
    """
    Creates an abject storing the properties of the Ionisation Chamber
    as defined in the settings file

    Parameters
    ----------
    settings : configparser
    Returns
    -------
    target : srim.target
    """
    layers = []
    preffix = "" #define a prefix to switch btw the simuation with and without the target
    if settings.getboolean("simulation", "target_in"):
        preffix = "with_target_"

    for i in range(settings.getint("target", preffix + "no_of_target_layers")): #loop over target layers
        elements_dict = {} #dictionary to store the elements of a single layer

        #define the composition of a layer
        for j in range(settings.getint("target", preffix + "layer_" + str(i) + "_no_of_elements")):
            element = settings.get("target", preffix + "layer_" + str(i) + "_element_" + str(j))
            element_dict = {"stoich" : settings.getint("target", preffix + "layer_" + str(i) + "_element_" + str(j) + "_stoich"),
                            "E_d" : settings.getfloat("target", element + "_displacement_energy"), # these values specify the energy losses characteristic for each element. Taken from the standalone SRIM 
                            "lattice" : settings.getfloat("target", element + "_lattice"),
                            "surface" : settings.getfloat("target", element + "_surface_energy")}
            elements_dict[element] = element_dict

        #check if the layer is gasous and define it appropriatelly
        #pysrim does not currently (27.08.25) support the compound correction directly. To circumvent, the density is scaled by the factor obtained from the standalone SRIM application. This overestimats the correction since it is equivalent to lowering of both the electronic and the nuclear stopping powers (only the former should be affected). This will introduce a small error (a fraction of a percent), dependending on the importance of nuclear interactions (which, for a high energy beam should be minimal).
        if settings.getboolean("target", preffix + "layer_" + str(i) + "_gas"):
            layers.append(Layer(elements_dict, density = gas_density_calc(settings.getfloat("target", preffix + "layer_" + str(i) + "_molar_mass") *
                                                                          settings.getfloat("target", preffix + "layer_" + str(i) + "_compound_correction"),
                                                                          settings.getfloat("target", preffix + "layer_" + str(i) + "_pressure"),
                                                                          temp = settings.getfloat("target", "temperature"),
                                                                          scaling_factor = 1 - simulation_uncertainty + settings.getfloat("simulation", "density_scaling_factor")), #density scaling factor and the simulation uncertaities go here
                                width = settings.getfloat("target", preffix + "layer_" + str(i) + "_thickness")))

        else:
            layers.append(Layer(elements_dict, density = settings.getfloat("target", preffix + "layer_" + str(i) + "_density") * (
                                simulation_uncertainty * settings.getfloat("target", preffix + "layer_" + str(i) + "_compound_correction") ), #simulation unceratinty added as a a scaling factor
                                width = settings.getfloat("target", preffix + "layer_" + str(i) + "_thickness")))
    return Target(layers) #define the target with the layers found

def define_beam(settings, iterator):
    """
    Creates an abject storing the properties of the target as defined in the settings file

    Parameters
    ----------
    settings : configparser
    iterator : int
        The index of the beam from the config file to be simulated 
    Returns
    -------
    ion : srim.ion

    """
    Z = settings.getint("beam", "beam_Z_" + str(iterator + 1))
    A = settings.getint("beam", "beam_A_" + str(iterator + 1))
    return Ion(elements[Z].symbol, energy = A * settings.getfloat("beam", "energy") * 10**6,
               mass = (A * atomic_mass_unit + mass.massExcess(Z, A)[0]) / atomic_mass_unit)


def plot_ionization(results, ax, ax_2, settings, atomic_number, mass_number):
    """
    Plots the energies deposited in the IC as a histogram

    Parameters
    ----------
    results : pandas dataframe
        Contains the raw results of the srim simulation
    ax : axes to plot on
    settings : configparser
    atomic_number : int
    mass_number : int
    Returns
    -------
    None.

    """
    #define the arrays for storing the results
    bins = []
    values = [] #for storing the integrated values
    values_binned = [[] for _ in range(13)] #integrated values for each segment to plot the projections
    E_signals = []
    dE_signals = []

    #loop over the simulated ions
    for i in range(settings.getint("simulation", "no_of_ions")):
        #variables for storing the E-dE signals
        dE_energy_start = 0
        dE_energy_stop = 0
        E_energy_start = 0
        E_energy_stop = 0

        ion_results = results[results["Ion#"] == i + 1] # extract results for a particular ion
        if len(ion_results) > 0 :
            #find the the distance to the first anode in Angstrom (add the target thickness as needed)
            start_pos = settings.getfloat("simulation", "dist_to_window") + settings.getfloat("simulation", "dist_to_anode") + (
                settings.getfloat("target", "with_target_layer_0_thickness") if settings.getboolean("simulation", "target_in") else 0 ) 
            #find index corresponding to the first entry in the region of interest
            start_index = (ion_results["Depth"] - start_pos).abs().idxmin()
            #loop over the anodes and get the energy that the ion has left in at the boundaries of the regions
            for j in range(13):
                stop_pos = start_pos + (j + 1) * settings.getfloat("simulation", "anode_spacing")
                stop_index = (ion_results["Depth"] - stop_pos).abs().idxmin() #index in the ion_results subset corresponding to the last entry in the region of interest
                #update results
                bins.append(j)
                values.append((ion_results["Energy"][start_index] - ion_results["Energy"][stop_index]) / 10**3)
                values_binned[j].append((ion_results["Energy"][start_index] - ion_results["Energy"][stop_index]) / 10**3)

                #get the E-dE values
                if j == settings.getint("plotting", "dE_start") :
                    dE_energy_start = ion_results["Energy"][start_index]
                elif j == settings.getint("plotting", "dE_stop") :
                    dE_energy_stop = ion_results["Energy"][stop_index]
                elif j == settings.getint("plotting", "E_start") :
                    E_energy_start = ion_results["Energy"][start_index]
                elif j == settings.getint("plotting", "E_stop") :
                    E_energy_stop = ion_results["Energy"][stop_index]

                start_index = stop_index #prepare for the next anode

        E_signals.append((E_energy_start - E_energy_stop) / 10**3)
        dE_signals.append((dE_energy_start - dE_energy_stop) / 10**3)

    #plot the results
    #The 2D hit Bragg curve
    ax[0].hist2d(bins, values, range = [[0, 13], [settings.getfloat("plotting", "bragg_energy_min"), settings.getfloat("plotting", "bragg_energy_max")]],
                 bins = [13, settings.getint("plotting", "bragg_energy_bins")],  cmap = "viridis", cmin = 1)
    ax[0].set_xlabel("Anode#")
    ax[0].set_ylabel("Energy (MeV)")
    #The dE-E hist
    ax[1].hist2d(dE_signals, E_signals, range = [[settings.getfloat("plotting", "dE_min"), settings.getfloat("plotting", "dE_max")],
                                                 [settings.getfloat("plotting", "E_min"), settings.getfloat("plotting", "E_max")]],
                 bins = [settings.getint("plotting", "dE_bins"), settings.getint("plotting", "E_bins")], cmap = "viridis", cmin = 1)
    ax[1].set_xlabel("dE (MeV)")
    ax[1].set_ylabel("E (MeV)")

    #the projections for individual segments
    for i in range(13):
        if i == 0:
            ax_2[i//4, i%4].hist(values_binned[i], bins=500, range = [0 , settings.getfloat("plotting", "projection_max")], label = str(
                atomic_number) + str(elements[atomic_number].symbol) + str(mass_number))
        else:
            ax_2[i//4, i%4].hist(values_binned[i], bins=500, range = [0 , settings.getfloat("plotting", "projection_max")])
        ax_2[i//4, i%4].set_xlabel("Energy (keV)")
        ax_2[i//4, i%4].set_ylabel("Counts")
        ax_2[i//4, i%4].set_title(f"Bragg curve from segment {i}")

    #save the data into txt .files to make it easier to compare with data
    data = np.column_stack((bins, values))
    filename = str(atomic_number) + str(elements[atomic_number].symbol) + str(mass_number) + ".txt"
    np.savetxt(filename, data, fmt="%.6f", header="Anode# Energy deposited")

    return 0

def error_analysis(results, settings):
    """
    Returns an array of average energies and their stds for each segment.
    Used for error estimation.

    Parameters
    ----------
    results : pandas dataframe
        Contains the results of the srim simulation
    settings : configparser
    
    Returns
    -------
    An array of average energies and their stds for each segment

    """

    #define the arrays for storing the results
    values_binned = [[] for _ in range(13)]

    #loop over the simulated ions
    for i in range(settings.getint("simulation", "no_of_ions")):

        ion_results = results[results["Ion#"] == i + 1] # extract results for a particular ion
        if len(ion_results) > 0 :
            #the distance to the first anode in Angstrom
            start_pos = settings.getfloat("simulation", "dist_to_window") + settings.getfloat("simulation", "dist_to_anode") + (
                settings.getfloat("target", "with_target_layer_0_thickness") if settings.getboolean("simulation", "target_in") else 0 )
            start_index = (ion_results["Depth"] - start_pos).abs().idxmin() #index in the ion_results subset corresponding to the first entry in the region of interest
            #loop over the anodes
            for j in range(13):
                stop_pos = start_pos + (j + 1) * settings.getfloat("simulation", "anode_spacing")
                stop_index = (ion_results["Depth"] - stop_pos).abs().idxmin() #index in the ion_results subset corresponding to the last entry in the region of interest

                values_binned[j].append((ion_results["Energy"][start_index] - ion_results["Energy"][stop_index]))

                start_index = stop_index #prepare for the next anode

    # Find the average and std for each segment
    centroids = []
    for i in range(len(values_binned)):
        #if np.average(values_binned[i]) < settings.getfloat("simulation", "signal_threshold"): #check if there is a signal is likely to be picked up by the detector (signal threshold)
        #    centroids.append(0)
        #    centroids.append(0)
        #else:
        centroids.append(np.average(values_binned[i]))
        centroids.append(np.std(values_binned[i]))

    return centroids

def get_energy_resolution(energy):
    """
    Returns the predicted energy resolution in keV (as an std of a gaussian fit) for
    a given energy using an empirical data obtained during the 2025 experimental campaign
    May slighty overestimate the std for an optimal situation (some mistakes during setup of the detector were made)
    Formula tested only until about 120MeV

    Parameters
    ----------
    energy : float in keV
    ----------
    Returns
    energy resolution in keV
    """
    return 500 + 0.016 * energy

def get_energy_resolution_at_last_segment(energy):
    """
    Returns the predicted energy resolution in keV at the last excited segment of the IC (as an std of a gaussian fit) for
    a given average energy measured at the last two semgent using an empirical data obtained during the 2025 experimental campaign
    The error on the fit is considerable (uncertainty of the order of 8% should be expected).
    May slighty overestimate the std for an optimal situation (some mistakes during setup of the detector were made)
    Formula tested only until about 50MeV

    Parameters
    ----------
    energy : float average energy deposited at the last two segments of the IC in keV
    ----------
    Returns
    energy resolution in keV
    """
    if energy < 20000: #to avoid problems with a tiny energy deposited in the last segment, which may lead to an unrealistically low predicted resolution
        return 1800
    return 1170 + 0.031 * energy

def compute_overlap(beams_centroids, pair_to_check, ax, settings):
    """
    Finds the percentage overlap between the two gaussian distributed spectra
    with given means (and stds taken from a empirical formula)

    Parameters
    ----------
    beams_centroids : [float]
        The centroids and the stds for each segment of each beams 
    fraction_1, fraction_2 : float 
        The fraction of the beam corresponding to each of the two components
    ax : axes to plot on (for visualisation purposes)
    """
    
    #find the number of segments with some energy deposited for each beam (ignore the segments with a tiny signal, this is a somewhat arbitrary number, but something in that ballpark is needed to sensibly differentiate between the two resolution regimes sensibly)
    range_idx_1 = np.count_nonzero(np.where(np.array(beams_centroids[2 * pair_to_check[0]][2:-1:2]) < settings.getfloat("simulation", "signal_threshold"), 0 , np.array(beams_centroids[2 * pair_to_check[0]][2:-1:2]))) - 1 # the indices at which the last bit of energy is deposited
    range_idx_2 = np.count_nonzero(np.where(np.array(beams_centroids[2 * pair_to_check[1]][2:-1:2]) < settings.getfloat("simulation", "signal_threshold"), 0 , np.array(beams_centroids[2 * pair_to_check[1]][2:-1:2]))) - 1
    #compute the maximum separation in ennergy between the two beams in a corresponding segment
    absolute_separation = abs(np.array(beams_centroids[2 * pair_to_check[0]][2:-1:2]) - np.array(beams_centroids[2 * pair_to_check[1]][2:-1:2]))
    #find the resolution of each beam at the each segment
    mean_energy_1 = np.average(beams_centroids[2 * pair_to_check[0]][2:-4:2]) # the average energy across the first few segments used to predict the resolution
    mean_energy_2 = np.average(beams_centroids[2 * pair_to_check[1]][2:-4:2])
    resolutions_1 = np.array([])
    resolutions_2 = np.array([])
    
    for i in range(len(absolute_separation)):
        if i >= range_idx_1 - 1: # if the segment with maximum separation is comes from either of the last two segments the resolution should be estimated using the formula accounting for energy straggling
            mean_energy_1 = np.average(beams_centroids[2 * pair_to_check[0]][2:-1:2][range_idx_1 - 1 : range_idx_1 + 1])
            resolutions_1 = np.insert(resolutions_1, len(resolutions_1), get_energy_resolution_at_last_segment(mean_energy_1))
            print(beams_centroids[2 * pair_to_check[0]][2:-1:2][i], resolutions_1[-1])
        else:
            mean_energy_1 = np.average(beams_centroids[2 * pair_to_check[0]][2:-4:2])
            resolutions_1 = np.insert(resolutions_1, len(resolutions_1), get_energy_resolution(mean_energy_1))
            print(beams_centroids[2 * pair_to_check[0]][2:-1:2][i], resolutions_1[-1])
        if i >= range_idx_2 - 1:
            mean_energy_2 = np.average(beams_centroids[2 * pair_to_check[1]][2:-1:2][range_idx_2 - 1 : range_idx_2 + 1])
            resolutions_2 = np.insert(resolutions_2, len(resolutions_2), get_energy_resolution_at_last_segment(mean_energy_2))
            print(beams_centroids[2 * pair_to_check[1]][2:-1:2][i], resolutions_2[-1])
        else:
            mean_energy_2 = np.average(beams_centroids[2 * pair_to_check[1]][2:-4:2])
            resolutions_2 = np.insert(resolutions_2, len(resolutions_2), get_energy_resolution(mean_energy_2))
            print(beams_centroids[2 * pair_to_check[1]][2:-1:2][i], resolutions_2[-1])
        
    #figure out the how well separated the beams are given their predicted resolution
    max_separation_idx = np.argmax(absolute_separation / (np.array(resolutions_1) + np.array(resolutions_2) ))
    
    #get the relevant information about the beams and their predicted resolution
    fraction_1 = settings.getfloat("beam", "beam_fraction_" + str(pair_to_check[0] + 1))
    fraction_2 = settings.getfloat("beam", "beam_fraction_" + str(pair_to_check[1] + 1))
    
    #mean energies of the peaks at the segment where the separation between the beams of interest is the largest
    energy_1 = beams_centroids[2 * pair_to_check[0]][2:-1:2][max_separation_idx]
    energy_2 = beams_centroids[2 * pair_to_check[1]][2:-1:2][max_separation_idx]

    #define functions used for plotting an computing the overlap
    def gaussian(x, energy, resolution, fraction):
        return fraction * np.exp(-(x - energy)**2 / (2 * resolution ** 2)) / (np.sqrt(2 * np.pi) * resolution)
    def gaussian_min(x):
        return np.minimum(gaussian(x, energy_1, resolutions_1[max_separation_idx], 1), gaussian(x, energy_2, resolutions_2[max_separation_idx], fraction_2 / (fraction_1 + fraction_2) ))
    
    result, err = quad(gaussian_min, min(energy_1 - 5 * resolutions_1[max_separation_idx], energy_2 - 5 * resolutions_2[max_separation_idx]),
                       max(energy_1 + 5 * resolutions_1[max_separation_idx], energy_2 + 5 * resolutions_2[max_separation_idx]))
    
    print("The maximum energy separation of", elements[beams_centroids[2 * pair_to_check[0] + 1][0]].symbol, beams_centroids[2 * pair_to_check[0] + 1][1], 
              "and", elements[beams_centroids[2 * pair_to_check[1] + 1][0]].symbol, beams_centroids[2 * pair_to_check[1] + 1][1], "beams is estimated to be",
              "{:.2f}".format(abs(energy_1 - energy_2)), "keV at segment", max_separation_idx + 1, ". With resolutions of", "{:.2f}".format(resolutions_1[max_separation_idx]), "and", "{:.2f}".format(resolutions_2[max_separation_idx]),
              "keV, respectively, the overlap between their spectra is estimated to be", "{:.2f}".format(result * 100), "%." )
    
    
    x = np.linspace(min(energy_1 - 5 * resolutions_1[max_separation_idx], energy_2 - 5 * resolutions_2[max_separation_idx]),
                    max(energy_1 + 5 * resolutions_1[max_separation_idx], energy_2 + 5 * resolutions_2[max_separation_idx]), 1000)
    
    label_1 = elements[beams_centroids[2 * pair_to_check[0] + 1][0]].symbol + " " + str(beams_centroids[2 * pair_to_check[0] + 1][1])
    label_2 = elements[beams_centroids[2 * pair_to_check[1] + 1][0]].symbol + " " + str(beams_centroids[2 * pair_to_check[1] + 1][1])
    
    ax.plot(x, gaussian(x, energy_1, resolutions_1[max_separation_idx], fraction_1), label = label_1, color = "r")
    ax.plot(x, gaussian(x, energy_2, resolutions_2[max_separation_idx], fraction_2), label = label_2, color = "g")
    ax.plot(x, gaussian(x, energy_1, resolutions_1[max_separation_idx], fraction_1) + gaussian(x, energy_2, resolutions_2[max_separation_idx], fraction_2), label = "Combined", color = "b")
    ax.legend()
    
    return max_separation_idx

def main():
    """
    Main function performing the simulation and plotting the results
    """
    os.chdir(os.path.dirname(os.path.abspath(__file__))) #change the working directory to the location of the script to make sure the config file is found (the code should work without this in most cases)
    settings = get_settings("config.txt") #get the settings
    targets = []
    if settings.getboolean("simulation", "error_analysis"): #allow for some density variations of the targets to simulate the uncertainty interval
        targets.append(define_target(settings)) #the nominal target (no simulation uncertanties)
        uncertainty = settings.getfloat("simulation", "simulation_uncertainty")
        for i in range( int( (uncertainty - 1) // 0.01 ) ): #introduce the uncertainty in steps of 0.01 on the density correction factor (when just edge values are taken, it may underestimate the uncertaint in the segment where the three curves meet)
            targets.append(define_target(settings, simulation_uncertainty = uncertainty - i * 0.01))
            targets.append(define_target(settings, simulation_uncertainty = 1 - i * 0.01))

    else:
        targets.append(define_target(settings)) #setup the target

    # initialise the dE-E plot
    fig, ax = plt.subplots(2,1, figsize = (9, 10))
    beam_name = ""
    for i in range(settings.getint("beam", "beam_components")):
        beam_name += str(settings.getint("beam", "beam_Z_" + str(i + 1))) + str(elements[settings.getint("beam", "beam_Z_" + str(i + 1))].symbol) + str(settings.getint("beam", "beam_A_" + str(i + 1)))
        if i != settings.getint("beam", "beam_components") - 1:
            beam_name += ", "
    ax[0].set_title("The Bragg curve in the IC for the " + beam_name + " beam (uncertainties underestimated)")
    ax[1].set_title("The dE-E spectrum for the " + beam_name + " beam (uncertainties underestimated)")
    # initialise the spectra plot
    fig_2, ax_2 = plt.subplots(4, 4, figsize = (15, 15), constrained_layout = True)
    #initialise the error plot
    fig_3, ax_3 = plt.subplots(1, 1, figsize = (10, 6))
    if not settings.getboolean("simulation", "error_analysis"): #check if necessary
        plt.close(fig_3)

    #loop across all of the relevant beams
    beams_centroids = []
    for i in range(settings.getint("beam", "beam_components")):
        error_centroids = [] #initialise for storing the located centroids
        for j in range(len(targets)): #loop over the targets to get the uncertainty
            if (settings.getboolean("simulation", "On") or len(targets) > 1) :
                beam = define_beam(settings, i) #setup the beam
                # Initialize a TRIM calculation
                trim = TRIM(targets[j], beam, number_ions=settings.getint("simulation", "no_of_ions"), calculation=1, collisions = True)
                srim_executable_directory = settings.get("simulation", "SRIM_directory") # Specify the directory of SRIM.exe
                trim.run(srim_executable_directory) #run the simulation

            # Load the results
            results = pd.read_csv(settings.get("simulation", "SRIM_directory") + "/SRIM Outputs/COLLISON.txt" , sep=r"³",
                                 skiprows=(76 if settings.getboolean("simulation", "target_in") else 72), header=None, encoding="latin1", on_bad_lines='skip',
                                 engine = 'python', skipinitialspace = True)#the number of lines to skip depends on the number of target layers

            results = results.drop(results.columns[[0, 4, 5, 6, 7, 8, 9, 10, 11, 12]], axis=1) #drop useless columns
            results.columns = ["Ion#", "Energy", "Depth"] # Give meaningful column names

            if j == 0:
                #plot the results for the 1st nominal target (no errors on the target)
                plot_ionization(results, ax, ax_2, settings, settings.getint("beam", "beam_Z_" + str(i + 1)), settings.getint("beam", "beam_A_" + str(i + 1)))
                beams_centroids.append(error_analysis(results, settings))
                beams_centroids.append( ( settings.getint("beam", "beam_Z_" + str(i + 1)) , settings.getint("beam", "beam_A_" + str(i + 1)) ) )

            if len(targets) > 1: #Find the centroids foe the error analysis
                error_centroids.append(error_analysis(results, settings))

        if len(error_centroids) > 1: 
            # Find the maximum and minimum value at each segment
            top_boundary = [] #for storing the limits
            bot_boundary = []
            optimal_values = []
            bins = []
            for j in range(int(len(error_centroids[0]) / 2)):
                #locate the boundares
                values = []
                for k in range(len(error_centroids)):
                    values.append(error_centroids[k][2 * j])
                bins.append(j)
                top_boundary.append(max(values))
                if max(values) != min(values):
                    bot_boundary.append(min(values))
                else:
                    bot_boundary.append(0)

                optimal_values.append(values[0])

            #plot the uncertainty intervals and add legend
            ax_3.plot(bins, optimal_values, label = str(settings.getint("beam", "beam_Z_" + str(i + 1))) + str(elements[settings.getint("beam", "beam_Z_" + str(i + 1))].symbol) + str(settings.getint("beam", "beam_A_" + str(i + 1))),
                      marker = 'o', markersize = 4,
                linestyle = '-', linewidth = 1, color = ['0','g','r','v'][i % 4])
            #labelling of the graph
            ax_3.set_xlabel("Segment index")
            ax_3.set_ylabel("Energy deposited per segment (KeV)")
            ax_3.set_title("The uncertainty intervals for energy deposited in each segment of the IC")
            ax_3.legend()

            ax_3.fill_between(bins, top_boundary, bot_boundary, alpha = 0.4)

            #save the error data into .txt file
            data_to_export = np.column_stack((bins, bot_boundary, top_boundary))
            filename = str(settings.getint("beam", "beam_Z_" + str(i + 1))) + str(elements[settings.getint("beam", "beam_Z_" + str(i + 1))].symbol) + str(settings.getint("beam", "beam_A_" + str(i + 1))) + "_errors.txt"
            np.savetxt(filename, data_to_export, fmt="%.6f", header="Anode# Min Max")

    pairs_to_check = [] #find unique beam combinations to figure out if they are resolvable
    for i in range(int(len(beams_centroids) / 2)):
        for j in range(int(len(beams_centroids) / 2)):
            if (sorted( [i,j] ) not in pairs_to_check and i != j):
                pairs_to_check.append(sorted( [i,j] ))
    
    axs = []
    for i in range(len(pairs_to_check)): #for each pair of beams, figure out if they are resolvable using the empirical energy resolution and compute the overlap
        fig, ax = plt.subplots(1, 1, figsize = (10, 6))
        axs.append(ax)
        compute_overlap(beams_centroids, pairs_to_check[i], ax, settings)
        
        # settings for the plots
        ax.set_xlabel("Energy deposited (keV)")
        ax.set_ylabel("Relative counts")
        ax.set_title("The overlap of the " + elements[beams_centroids[2 * pairs_to_check[i][0] + 1][0]].symbol + " " + str(beams_centroids[2 * pairs_to_check[i][0] + 1][1]) +
                        " and " + elements[beams_centroids[2 * pairs_to_check[i][1] + 1][0]].symbol + " " + str(beams_centroids[2 * pairs_to_check[i][1] + 1][1]) +
                        " beams at the segment with the maximum energy separation")
        
    #display the results
    fig_2.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

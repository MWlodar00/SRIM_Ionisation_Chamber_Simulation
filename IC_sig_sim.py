# -*- coding: utf-8 -*-
"""
Michał Włodarczyk

Script automating simulations of the signals for the Ionisation Chamber with SRIM.

"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from srim import TRIM, Ion, Layer, Target
from nuclyr import mass
from periodictable import elements

import configparser
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

def gas_density_calc(molar_mass, pressure, temp = 293) :
    """
    Get the density of a gas given its molar mass and pressure, assuming
    room temperature and using the ideal gas law

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
    return molar_mass * pressure / (temp * gas_const) / 10**4

def define_target(settings):
    """
    Creates an abject storing the properties of the Ionisation Chamber
    as defined in the settings file

    Parameters SB 15 09
    ----------
    settings : configparser
    Returns
    -------
    target : srim.target
    """
    layers = []
    for i in range(settings.getint("target", "no_of_target_layers")):#loop over target layers
        elements_dict = {} #dictionary to store the elements of a single layer

        #define the composition of a layer
        for j in range(settings.getint("target", "layer_" + str(i) + "_no_of_elements")):
            element = settings.get("target", "layer_" + str(i) + "_element_" + str(j))
            element_dict = {"stoich" : settings.getint("target", "layer_" + str(i) + "_element_" + str(j) + "_stoich"),
                            "E_d" : settings.getfloat("target", element + "_displacement_energy"), # these values specify the energy losses characteristic for each element. They were obtained from the standalone SRIM (this version uses defaults that produce incorrect results)
                            "lattice" : settings.getfloat("target", element + "_lattice"),
                            "surface" : settings.getfloat("target", element + "_surface_energy")}
            elements_dict[element] = element_dict

        #check if the layer is gasous and define it appropriatelly
        #pysrim does not currently (27.08.25) support the compound correction directly. To circumvent, the density is scaled by the factor obtained from the standalone SRIM application. This overestimats the correction since it is equivalent to lowering of both the electronic and the nuclear stopping powers (only the former should be affected). This will introduce a small error (a fraction of a percent), dependending on the importance of nuclear interactions (which, for a high energy beam should not be significant). The proper fix would require modifying the pysrim files (the layer and compound classes).
        if settings.getboolean("target", "layer_" + str(i) + "_gas"):
            layers.append(Layer(elements_dict, density = gas_density_calc(settings.getfloat("target", "layer_" + str(i) + "_molar_mass") *
                                                                          settings.getfloat("target", "layer_" + str(i) + "_compound_correction"),
                                                                          settings.getfloat("target", "layer_" + str(i) + "_pressure"),
                                                                          temp = settings.getfloat("target", "temperature")),
                                width = settings.getfloat("target", "layer_" + str(i) + "_thickness")))

        else:
            layers.append(Layer(elements_dict, density = settings.getfloat("target", "layer_" + str(i) + "_density") *
                                settings.getfloat("target", "layer_" + str(i) + "_compound_correction"),
                                width = settings.getfloat("target", "layer_" + str(i) + "_thickness")))
    return Target(layers) #define the target with the layers found

def define_beam(settings, iterator):
    """
    Creates an abject storing the properties of the target as defined in the settings file

    Parameters
    ----------
    settings : configparser
    Z_offset : int
        The atomic number beam will be this many units away from the one given in the settings file
    A_offset : int
        The mass number beam will be this many units away from the one given in the settings file
    Returns
    -------
    ion : srim.ion

    """
    Z = settings.getint("beam", "beam_Z_" + str(iterator + 1))
    A = settings.getint("beam", "beam_A_" + str(iterator + 1))
    return Ion(elements[Z].symbol, energy = A * settings.getfloat("beam", "energy") * 10**6,
               mass = (A * atomic_mass_unit + mass.massExcess(Z, A)[0]) / atomic_mass_unit)


def plot_ionization(results, ax, settings, atomic_number, mass_number):
    """
    Plots the energies deposited in the IC as a histogram

    Parameters
    ----------
    results : pandas dataframe
        Contains the results of the srim simulation
    ax : asix
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
    E_signals = []
    dE_signals = []

    #loop over the simulated ions
    for i in range(settings.getint("simulation", "no_of_ions")):
        #variables for storing the energies deposited at particular anodes (for E-dE signals)
        first_anode_energy = 0
        last_anode_energy = 0
        intermediate_anode_energy = 0

        ion_results = results[results["Ion#"] == i + 1] # extract results for a particular ion
        if len(ion_results) > 0 :
            start_pos = settings.getfloat("simulation", "dist_to_window") + settings.getfloat("simulation", "dist_to_anode")  #the distance to the first anode in Angstrom
            start_index = (ion_results["Depth"] - start_pos).abs().idxmin() #index in the ion_results subset corresponding to the first entry in the region of interest
            #loop over the anodes
            for j in range(13):
                stop_pos = start_pos + (j + 1) * settings.getfloat("simulation", "anode_spacing")
                stop_index = (ion_results["Depth"] - stop_pos).abs().idxmin() #index in the ion_results subset corresponding to the last entry in the region of interest
                #update results
                bins.append(j)
                values.append((ion_results["Energy"][start_index] - ion_results["Energy"][stop_index]) / 10**3)

                #get the E-dE values
                if j == 0 :
                    first_anode_energy = ion_results["Energy"][start_index]
                elif j == 13 :
                    last_anode_energy = ion_results["Energy"][stop_index]
                elif j == settings.getfloat("simulation", "intermediate_anode_index") - 1 :
                    intermediate_anode_energy = ion_results["Energy"][start_index]

                start_index = stop_index #prepare for the next anode

        E_signals.append((first_anode_energy - last_anode_energy) / 10**3)
        dE_signals.append((intermediate_anode_energy - last_anode_energy) / 10**3)

    #plot the results
    ax[0].hist2d(bins, values, range = [[0, 13], [settings.getfloat("plotting", "bragg_energy_min"), settings.getfloat("plotting", "bragg_energy_max")]],
                 bins = [13, settings.getint("plotting", "bragg_energy_bins")],  cmap = "viridis", cmin = 1)
    ax[0].set_xlabel("Anode#")
    ax[0].set_ylabel("Energy (MeV)")
    ax[1].hist2d(dE_signals, E_signals, range = [[settings.getfloat("plotting", "dE_min"), settings.getfloat("plotting", "dE_max")],
                                                 [settings.getfloat("plotting", "E_min"), settings.getfloat("plotting", "E_max")]],
                 bins = [settings.getint("plotting", "dE_bins"), settings.getint("plotting", "E_bins")], cmap = "viridis", cmin = 1)
    ax[1].set_xlabel("dE (MeV)")
    ax[1].set_ylabel("E (MeV)")

    #save the data into txt .files
    data = np.column_stack((bins, values))
    filename = str(atomic_number) + str(elements[atomic_number].symbol) + str(mass_number) + ".txt"
    np.savetxt(filename, data, fmt="%.6f", header="Anode# Energy deposited")

def main():
    """
    Main function performing the simulation and plotting the results
    """
    settings = get_settings("config.txt") #get the settings
    target = define_target(settings) #setup the target

    # initialise the plot
    fig, ax = plt.subplots(2,1, figsize = (6, 10))

    #prepare the name of the plot
    name = str(settings.getint("beam", "beam_Z_1")) + str(elements[settings.getint("beam", "beam_Z_1")].symbol) + str(
        settings.getint("beam", "beam_A_1"))

    #loop across all of the relevant beams
    for i in range(settings.getint("beam", "beam_components")):
        if (settings.getboolean("simulation", "On")) :
            beam = define_beam(settings, i) #setup the beam
            # Initialize a TRIM calculation
            trim = TRIM(target, beam, number_ions=settings.getint("simulation", "no_of_ions"), calculation=1, collisions = True)
            srim_executable_directory = settings.get("simulation", "SRIM_directory") # Specify the directory of SRIM.exe
            trim.run(srim_executable_directory) #run the simulation

        # Load the results
        results = pd.read_csv(settings.get("simulation", "SRIM_directory") + "/SRIM Outputs/COLLISON.txt" , sep=r"³",
                             skiprows=72, header=None, encoding="latin1", on_bad_lines='skip',
                             engine = 'python', skipinitialspace = True)

        results = results.drop(results.columns[[0, 4, 5, 6, 7, 8, 9, 10, 11, 12]], axis=1) #drop useless columns
        results.columns = ["Ion#", "Energy", "Depth"] # Give meaningful column names
        #plot the results
        plot_ionization(results, ax, settings, settings.getint("beam", "beam_Z_" + str(i + 1)), settings.getint("beam", "beam_A_" + str(i + 1)))
        #update the name if needed
        if (i != 0) :
            name += ", " + str(settings.getint("beam", "beam_Z_" + str(i + 1))) + elements[settings.getint(
                "beam", "beam_Z_"+ str(i + 1))].symbol + str(settings.getint("beam", "beam_A_" + str(i + 1)))

    #display the results
    ax[0].set_title(name)
    plt.tight_layout()

main()

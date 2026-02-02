# Data to simulation comparison
#NOTES: For gunwati's hdf -> .npy burnsample

# imports
import numpy as np
import matplotlib.pyplot as plt
import math
import h5py, glob
import sys, os, glob, argparse, ast
from ang_res_funcs import *
from itertools import chain


#import pandas as pt

def load_simulation(input_directory, param, tier, model):
    # Zenith
    particles = ['Proton', 'Helium', 'Oxygen', 'Iron']
    all_part = []
    weights = []
    for p in particles:
        if param == "Stations":
            all_part.append(np.load(input_directory+'/{}-{}-{}.npy'.format(p, tier, param)))
        elif param == "Energy":
            all_part.append(np.load(input_directory+'/{}-{}-MC-{}.npy'.format(p, tier, param)))
        else:
            all_part.append(np.load(input_directory+'/{}-{}-{}-{}.npy'.format(p, tier, model, param)))
        weights.append(np.load(input_directory+'/{}-{}-{}-Weights.npy'.format(p, tier, model)))
        
    output_one = np.concatenate((all_part[0], all_part[1], all_part[2], all_part[3]), axis=None)
    output_two = np.concatenate((weights[0], weights[1], weights[2], weights[3]), axis=None)

    return [output_one, output_two]
    
def load_burnsample(input_file):
    dat = np.load(input_file)
    stations = dat[1]
    zenith = dat[2]
    azimuth = dat[3]
    return [stations, zenith, azimuth]
        
def apply_cuts(arr, stations, bounds, radians=False):
    arrays = []
    for band in bounds:
        tier = arr[np.where((stations >= band[0]) & (stations < band[1]))]
        arrays.append(tier)
    final = arr[stations >= band[1]]
    arrays.append(final)

    output = []
    #Convert to degrees
    if radians:
        for tier in arrays:
            output.append(np.rad2deg(tier))
        return output
    else:
        return arrays

if __name__ == '__main__':
    # Set up command line options
    parser = argparse.ArgumentParser(description=
                                     "Compares data to simulation.")
    #Must supply first two arguments
    parser.add_argument('data_file', type=str,
                        help="Path of burnsample .npy\
                        GA default: '/data/user/gagrawal/arrays/it73/att2/burnsample_2012.npy'")
    parser.add_argument('sim_directory', type=str,
                        help="Directory where simulation .npy files are stored.\
                        GA default: '/home/gagrawal/stages/jsummers/arrays7'")
    # One of the two
    parser.add_argument("-z", "--plot_zenith", dest="plot_zenith", default=False, 
                        action="store_true",
                        help="Produces plots for zenith distribution.")
    parser.add_argument("-a", "--plot_azimuth", dest="plot_azimuth", default=False,
                        action="store_true",
                        help="Produces plots for azimuth distribution.")
    parser.add_argument("--stations", dest="stations", action="store_true", default=False,
                        help="NStation distribution.")
    parser.add_argument("--years", dest="stations_sim", type=str, nargs=3,
                        default='', #default='./arrays/default ./arrays/2015/ns ./arrays/2018/ns',
                        help="Scaled NStation distribution using 3 years of simulation.")
    #Optional
    parser.add_argument("-o", "--output", dest="output", type=str,
                        default=".",
                        help="Output directory of plots")
    parser.add_argument("-e", "--energies", dest='energies', type=str,
                        default='[[3, 5], [5, 9], [9, 16]]',
                        help="NStation boundaries for each energy tier/band for 2012")
    parser.add_argument("-n", "--number", dest='number', type=int,
                        default=15, help="Number of (linear) bins")
    parser.add_argument("-r", "--ratio", dest="ratio", action="store_true",
                        default=False, help="Plot ratio of sim/data.")
    parser.add_argument("--log", dest="log", action="store_true", default=False,
                        help="Log scale.")
    parser.add_argument("--cos", dest="cos_bin", action="store_true", default=False,
                        help="Cosine binning.")
    args = parser.parse_args()
        
    # Parse energies as list of lists
    if args.energies is not None:
        try:
            bins = ast.literal_eval(args.energies)
            if not all(isinstance(inner_list, list) and all(isinstance(i, int) for i in inner_list) for inner_list in bins):
                raise ValueError
        except (ValueError, SyntaxError):
            print(f"Error: The provided list of lists '{args.energies}' is not valid.")
            sys.exit(1)

    # if output directory does not exist, create it
    if not os.path.isdir(args.output):
        print('Creating output directory {}'.format(args.output))
        os.makedirs(args.output)

    reco_key = [(1, "ShowerPlane"), (2, "Laputop"), (3, "Laputop"), (4, "Laputop")]
    
    # Load data sets
    stations, zenith, azimuth = load_burnsample(args.data_file)

    if args.plot_zenith:
        if args.cos_bin:
            bin_ends = get_cosine_binning(62.859, 11)
        else:
            bin_ends = np.linspace(0, 60, args.number+1)
        # Bin burnsample with the energies
        data = []
        tiers = apply_cuts(zenith, stations, bins, True)
        for tier in tiers:
            dat_err = bin_weighted_param(tier, bin_ends)[-1]
            data.append([tier, dat_err])
        
        #Load simulations
        simulations = []
        for ind, reco in reco_key:
            sim_zenith, sim_weights = load_simulation(args.sim_directory, "Zenith", "T{}".format(ind), reco)
            err = bin_weighted_param(sim_zenith, bin_ends, sim_weights)[-1]
            simulations.append([sim_zenith, sim_weights, err])

        if args.ratio:
            ratio_sim_data_all_tiers(data, simulations, "Zenith", bin_ends,
                                     reco_key, logscale=args.log, output=args.output+'/ratio-2012-zenith.png')
        else:
            dist_sim_data_all_tiers(data, simulations, "Zenith", bin_ends,
                                    reco_key, logscale=args.log, output=args.output+'/simvdat-2012-zenith.png')

    if args.plot_azimuth:
        bin_ends = np.linspace(0, 360, args.number+1)
        # Bin burnsample with the energies
        data = []
        tiers = apply_cuts(azimuth, stations, bins, True)
        for tier in tiers:
            dat_err = bin_weighted_param(tier, bin_ends)[-1]
            data.append([tier, dat_err])
            
        #Load simulations
        simulations = []
        for ind, reco in reco_key:
            sim_azimuth, sim_weights = load_simulation(args.sim_directory, "Azimuth", "T{}".format(ind), reco)
            err = bin_weighted_param(sim_azimuth, bin_ends, sim_weights)[-1]
            simulations.append([sim_azimuth, sim_weights, err])

        if args.ratio:
            ratio_sim_data_all_tiers(data, simulations, "Azimuth", bin_ends,
                                     reco_key, logscale=args.log, output=args.output+'/ratio-2012-azimuth.png')
        else:
            dist_sim_data_all_tiers(data, simulations, "Azimuth", bin_ends,
                                    reco_key, logscale=args.log, output=args.output+'/simvdat-2012-azimuth.png')

    if args.stations:
        bottom = 5
        top = 22
        bin_ends = np.arange(bottom, top, 1)
        bin_centres = bin_ends[:-1]        

        sim_nstat = load_simulation(args.sim_directory, "Stations", "T1", "MC")[0]
        sim_weights = []
        tiers = apply_cuts(sim_nstat, sim_nstat, bins, radians=False)
        for ind, reco in reco_key:
            sim_weights.append(load_simulation(args.sim_directory, "Stations", "T{}".format(ind), "MC")[-1])

        sim_ns = np.array(list(chain.from_iterable(tiers)))
        sim_w = np.array(list(chain.from_iterable(sim_weights)))
        
        err = bin_weighted_param(sim_ns, bin_ends, sim_w)[-1]
        y_reco = np.histogram(sim_ns, bin_ends, weights=sim_w)[0]
        plt.errorbar(bin_centres, y_reco/sum(y_reco), yerr=err/sum(y_reco), label='Simulation', capsize=3)

        data, data_err = [stations, bin_weighted_param(stations, bin_ends)[-1]]
        y_data = np.histogram(data, bin_ends)[0]
        plt.errorbar(bin_centres, y_data/sum(y_data), yerr=data_err/sum(y_data), label='Burnsample', capsize=3)
        
        plt.ylabel('Fraction of Events', fontsize=15)
        plt.xlabel('NStations', fontsize=15)
        plt.xlim([bottom - 1, top - 1])
        plt.xticks(np.arange(bottom, top - 1, 1))
        plt.title('Sim vs Data NStation Distribution', fontsize=15)
        plt.legend(loc='upper right')
        
        plt.savefig('./sim_data_nS_2012.png', transparent=False)
   
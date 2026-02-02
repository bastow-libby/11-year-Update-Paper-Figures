#imports
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import pandas as pd
from itertools import product

import math
import h5py, glob
import sys, os, glob, argparse, ast
from ang_res_funcs import *
from itertools import chain

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

     
def apply_cuts(arr, stations, bounds):
    arrays = []
    for band in bounds:
        tier = arr[np.where((stations >= band[0]) & (stations < band[1]))]
        arrays.append(tier)
    final = arr[stations >= band[1]]
    arrays.append(final)

    return arrays

if __name__ == "__main__":
    reco_key = [(1, "ShowerPlane"), (2, "Laputop"), (3, "Laputop"), (4, "Laputop")]
    directories = ['./arrays/default', './arrays/2015/ns', './arrays/2018/ns']
    yearly = [2012, 2015, 2018]
    bins = [[3,5], [5,9], [9,16]]
    Colors = ['magenta', 'gold', 'cyan']
    bin_start = 5
    bin_end = 21
    
    bin_ends = np.arange(bin_start, bin_end+1, 1)
    bin_centres = bin_ends[:-1]
    fig = plt.figure(figsize=[6, 8])

    #For each simulation dataset
    for i, dirc in enumerate(directories):
        #Load all stations and split into tiers
        sim_nstat = load_simulation(dirc, "Stations", "T1", "MC")[0]
        tiers = apply_cuts(sim_nstat, sim_nstat, bins)

        #Load weights and zeniths for each tier
        sim_weights = []
        sim_zen = []
        for ind, reco in reco_key:
            sim_weights.append(load_simulation(dirc, "Weights", "T{}".format(ind), "MC")[0])
            sim_zen.append(load_simulation(dirc, "Zenith", "T{}".format(ind), "MC")[0])

        #Combine all tiers for each param
        sim_ns = np.array(list(chain.from_iterable(tiers)))
        sim_w = np.array(list(chain.from_iterable(sim_weights)))
        sim_z = np.array(list(chain.from_iterable(sim_zen)))

        #Zenith cut
        sim_ns = sim_ns[sim_z < 55]
        sim_w = sim_w[sim_z < 55]

        '''
        #weighted distribution 
        err = bin_weighted_param(sim_ns, bin_ends, sim_w)[-1]
        y_reco = np.histogram(sim_ns, bin_ends, weights=sim_w)[0]
        plt.errorbar(bin_centres, y_reco, yerr=err, label=yearly[i], capsize=3)
        '''
        
        #Bin
        y_reco = np.histogram(sim_ns, bin_ends, weights=sim_w)[0]
        counts = y_reco #/ sum(y_reco)
        if i==0:
            base = counts
            
        #plt.errorbar(bin_centres, counts/base, label=yearly[i], capsize=3)
        plt.stairs(counts/base, bin_ends, label=yearly[i], color=Colors[i], linewidth=2.0)
        
    plt.ylabel('Frac. of Counts / Frac. of Counts (2012)', fontsize=15)
    plt.xlabel('NStations', fontsize=15)
    plt.xlim(bin_start, bin_end)
    plt.xticks(np.arange(bin_start, bin_end, 1))
    plt.yticks(np.arange(0, 1.2, .1))
    plt.title('Sim NStation Distributions scaled to 2012', fontsize=15)
    plt.legend(loc='upper right')
        
    plt.savefig('./yelp.png', transparent=False)
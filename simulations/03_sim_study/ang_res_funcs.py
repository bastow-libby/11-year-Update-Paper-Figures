# imports
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pt

"""
This file contains methods that produce the angular resolution plots for our
analysis.


"""


# ANGULAR RESOLUTION
# Checking reconstruction accuracy by looking at difference between
# MCPrimary and Sim

def find_opening_angle(true_zen, reco_zen, true_az, reco_az):
    dlamda = np.deg2rad(reco_az - true_az)
    phi1 = np.deg2rad(90 - reco_zen)
    phi2 = np.deg2rad(90 - true_zen)
    return np.rad2deg(np.arccos(np.sin(phi1) * np.sin(phi2) +
                                np.cos(phi1) * np.cos(phi2) *
                                np.cos(dlamda)))


# We use cosine binning for the following plots

def get_cosine_binning(upper_bound, bin_count):
    angle_last = upper_bound * np.pi / 180

    cos_zenbin_ends = np.linspace(1, np.cos(angle_last), bin_count + 1)
    arccos_zenbin_ends = np.arccos(cos_zenbin_ends) * 180 / np.pi

    # Get center of bins
    arccos_zenbins = (arccos_zenbin_ends[1:] + arccos_zenbin_ends[:-1]) / 2

    # print(arccos_zenbin_ends)
    # print(arccos_zenbins)
    xerr = (arccos_zenbin_ends[1:] - arccos_zenbin_ends[:-1]) / 2

    # print(xerr)
    return arccos_zenbin_ends


def plot_ang_res_by_tier(data_list, bins_list, tiers, titles, param, x_coord):
    n = len(data_list)
    height = math.ceil(n / 2)
    fig = plt.figure(figsize=(12, 6 * height), constrained_layout=True)
    ax_array = fig.subplots(2, height, squeeze=False)
    fig.suptitle('Reconstruction Resolution for ' + param + ' By Tier',
                 fontsize=15)

    row = 0
    col = 0

    for data, bins, tier, title in zip(data_list, bins_list,
                                       tiers, titles):
        low_err = np.quantile(data, .16)
        high_err = np.quantile(data, .84)
        
        ax_array[row, col].hist(data, bins=bins)
        ax_array[row, col].set_title(tier, fontsize=15)
        ax_array[row, col].set_ylabel('Count', fontsize=15)
        ax_array[row, col].set_xlabel(title + ' ' + param, fontsize=15)
        ax_array[row, col].set_yscale('log')
        ax_array[row, col].vlines(np.mean(data), 0, 1, transform=ax_array[row, col].get_xaxis_transform(), linestyle='solid')
        ax_array[row, col].vlines(low_err, 0, 1, transform=ax_array[row, col].get_xaxis_transform(), linestyle='dashed')
        ax_array[row, col].vlines(high_err, 0, 1, transform=ax_array[row, col].get_xaxis_transform(), linestyle='dashed')
        
        s1 = 'Mean: {}'.format(round(np.mean(data), 3))
        s2 = 'RMS: {}'.format(round(np.sqrt(np.mean(data ** 2)), 3))
        s3 = '68%: {} - {}'.format(round(low_err, 3), round(high_err, 3))
        
        ax_array[row, col].annotate(s1, xy=(20, 60), xytext=(x_coord, 1000),
                                    fontsize=18)
        ax_array[row, col].annotate(s2, xy=(20, 40), xytext=(x_coord, 500),
                                    fontsize=18)
        ax_array[row, col].annotate(s3, xy=(20, 20), xytext=(x_coord, 200), fontsize=15)
        # ax_array[row, col].vlines(55, 0, 5000000, linestyle = 'dashed')

        col = col + 1
        if col == 2:
            col = 0
            row = row + 1


def reco_ang_res(true_zen, reco_zen, true_az, reco_az, reconstruction):
    # Plots histograms of zenith difference, azimuth difference, and opening
    # angle for Monte Carlo and a given reconstruction.

    zendiff = true_zen - reco_zen
    zbins = np.linspace(-60, 60, 121)
    zbincenters = 0.5 * (zbins[1:] + zbins[:-1])

    azdiff = true_az - reco_az
    abins = np.linspace(-360, 360, 721)
    abincenters = 0.5 * (abins[1:] + abins[:-1])

    opening_angle = find_opening_angle(true_zen, reco_zen, true_az, reco_az)

    fig = plt.figure(figsize=(15, 5), constrained_layout=True)
    ax_array = fig.subplots(1, 3, squeeze=False)
    fig.suptitle('Reconstruction Resolution for Zenith & Azimuth')

    ax_array[0, 0].hist(zendiff, zbincenters)
    ax_array[0, 0].set_title('Zenith')
    ax_array[0, 0].set_ylabel('Count')
    ax_array[0, 0].set_xlabel('MC Zenith - ' + reconstruction + ' Zenith')
    ax_array[0, 0].set_yscale('log')
    s1 = 'Mean: {}'.format(round(np.mean(zendiff), 3))
    s2 = 'RMS: {}'.format(round(np.sqrt(np.mean(zendiff ** 2)), 3))
    s3 = '68%: {} - {}'.format(round(np.quantile(zendiff, .16), 3), round(np.quantile(zendiff, .84), 3))
    ax_array[0, 0].annotate(s1, xy=(20, 60), xytext=(20, 1000))
    ax_array[0, 0].annotate(s2, xy=(20, 40), xytext=(20, 500))
    ax_array[0, 0].annotate(s3, xy=(20, 20), xytext=(23, 200))
    ax_array[0, 0].vlines(np.quantile(zendiff, .16), 0, 1, transform=ax_array[0, 0].get_xaxis_transform(), linestyle='dashed')
    ax_array[0, 0].vlines(np.mean(zendiff), 0, 1, transform=ax_array[0, 0].get_xaxis_transform(), linestyle='solid')
    ax_array[0, 0].vlines(np.quantile(zendiff, .84), 0, 1, transform=ax_array[0, 0].get_xaxis_transform(), linestyle='dashed')

    ax_array[0, 1].hist(azdiff, abincenters)
    ax_array[0, 1].set_title('Azimuth')
    ax_array[0, 1].set_ylabel('Count')
    ax_array[0, 1].set_xlabel('MC Azimuth - ' + reconstruction + ' Azimuth')
    ax_array[0, 1].set_yscale('log')
    s1 = 'Mean: {}'.format(round(np.mean(azdiff), 3))
    s2 = 'RMS: {}'.format(round(np.sqrt(np.mean(azdiff ** 2)), 3))
    s3 = '68%: {} - {}'.format(round(np.quantile(azdiff, .16), 3), round(np.quantile(azdiff, .84), 3))
    ax_array[0, 1].annotate(s1, xy=(10, 60), xytext=(100, 1000))
    ax_array[0, 1].annotate(s2, xy=(10, 40), xytext=(100, 500))
    ax_array[0, 1].annotate(s3, xy=(10, 20), xytext=(100, 200))
    ax_array[0, 1].vlines(np.quantile(azdiff, .16), 0, 1, transform=ax_array[0, 1].get_xaxis_transform(), linestyle='dashed')
    ax_array[0, 1].vlines(np.mean(azdiff), 0, 1, transform=ax_array[0, 1].get_xaxis_transform(), linestyle='solid')
    ax_array[0, 1].vlines(np.quantile(azdiff, .84), 0, 1, transform=ax_array[0, 1].get_xaxis_transform(), linestyle='dashed')


    ax_array[0, 2].hist(opening_angle, bins=100)
    ax_array[0, 2].set_title('Opening Angle')
    ax_array[0, 2].set_ylabel('Count')
    ax_array[0, 2].set_xlabel('Opening Angle')
    ax_array[0, 2].set_yscale('log')
    s1 = 'Mean: {}'.format(round(np.mean(opening_angle), 3))
    s2 = 'RMS: {}'.format(round(np.sqrt(np.mean(opening_angle ** 2)), 3))
    s3 = '68%: {} - {}'.format(round(np.quantile(opening_angle, .16), 3), round(np.quantile(opening_angle, .84), 3))
    ax_array[0, 2].annotate(s1, xy=(10, 60), xytext=(50, 1000))
    ax_array[0, 2].annotate(s2, xy=(10, 40), xytext=(50, 500))
    ax_array[0, 2].annotate(s3, xy=(10, 20), xytext=(50, 200))
    ax_array[0, 2].vlines(np.quantile(opening_angle, .16), 0, 1, transform=ax_array[0, 2].get_xaxis_transform(), linestyle='dashed')
    ax_array[0, 2].vlines(np.mean(opening_angle), 0, 1, transform=ax_array[0, 2].get_xaxis_transform(), linestyle='solid')
    ax_array[0, 2].vlines(np.quantile(opening_angle, .84), 0, 1, transform=ax_array[0, 2].get_xaxis_transform(), linestyle='dashed')


def reco_open_ang_res_bytier(mc_zen, reco_zen, mc_az, reco_az, reco_names, tiers):
    # Produce histograms of difference between MC and simulation in
    # Opening Angle for each energy tier.

    opening_angles = []
    obins = np.linspace(0, 180, 100)
    ob_centres = []

    for (true_zen, reco_zen, true_az, reco_az) in zip(mc_zen, reco_zen, mc_az,
                                                      reco_az):
        opening_angles.append(find_opening_angle(true_zen, reco_zen, true_az,
                                                 reco_az))
        ob_centres.append(0.5 * (obins[1:] + obins[:-1]))
    plot_ang_res_by_tier(opening_angles, ob_centres, tiers, reco_names,
                         "Opening Angle", 50)


def reco_zen_res_bytier(mc_zen, reco_zen, reco_names, tiers):
    zendiffs = []
    zbins = np.linspace(-60, 60, 121)
    titles = []
    zb_centres = []

    for (true_zen, reco_zen, reconstruction) in zip(mc_zen, reco_zen, reco_names):
        zendiffs.append(true_zen - reco_zen)
        zb_centres.append(0.5 * (zbins[1:] + zbins[:-1]))
        titles.append("MC Zenith - " + reconstruction)

    plot_ang_res_by_tier(zendiffs, zb_centres, tiers, titles, "Zenith", 20)


def reco_az_res_bytier(mc_az, reco_az, reco_names, tiers):
    azdiffs = []
    abins = np.linspace(-360, 360, 721)
    titles = []
    ab_centres = []

    for (true_az, reco_az, reconstruction) in zip(mc_az, reco_az, reco_names):
        azdiffs.append(true_az - reco_az)
        ab_centres.append(0.5 * (abins[1:] + abins[:-1]))
        titles.append("MC Azimuth - " + reconstruction)

    plot_ang_res_by_tier(azdiffs, ab_centres, tiers, titles, "Azimuth", 100)



# ZENITH DIFFERENCE PER ZENITH BIN
# Plot median zenith differnce between simulated and reconstructed events vs zenith angle

def plot_vs_zenith_bytier(data_list, titles, param):
    n = len(titles)
    height = math.ceil(n / 2)
    fig = plt.figure(figsize=(16, 4 * height), constrained_layout=True)
    ax_array = fig.subplots(2, height, squeeze=False)
    fig.suptitle('Median ' + param + ' Per Zenith Bin By Tier',
                 fontsize=15)

    row = 0
    col = 0

    for (data, title) in zip(data_list, titles):
        for particle in data:
            yave = particle[0]
            xerr = particle[1]
            yerr = particle[2]
            bins = particle[3]
            label = particle[4]
            color = particle[5]
            ax_array[row, col].errorbar(bins, yave, xerr=xerr, yerr=yerr, label=label, color=color, capsize=2,
                                        fmt='o')
        ax_array[row, col].vlines(55, 0, 2000, linestyle='dashed')
        ax_array[row, col].legend(loc='lower right')
        ax_array[row, col].set_title(title, fontsize=15)
        ax_array[row, col].set_ylabel(r'$\Delta\Theta\;[\degree]$', fontsize=12)
        ax_array[row, col].set_xlabel(r'$Zenith\;angle\;[\degree]$ ', fontsize=12)
        ax_array[row, col].set_ylim(.1, 100)
        ax_array[row, col].set_xlim(-2, 68)
        #Make tick marks larger!
        ax_array[row, col].set_xticks(np.arange(0, 68, 4))
        ax_array[row, col].set_yscale('log')

        col = col + 1
        if (col == 2):
            col = 0
            row = row + 1


def zenith_diff(true_zen, reco_zen, bin_ends, cut=65):
    # Zenith difference

    dzen = np.abs(reco_zen - true_zen)
    # The simulation does not include any zenith values above 65 degrees, so we will cut out values that get too high
    y = np.delete(dzen, np.argwhere(reco_zen >= cut))
    x = np.delete(reco_zen, np.argwhere(reco_zen >= cut))

    # Get center of bins
    bins = (bin_ends[1:] + bin_ends[:-1]) / 2
    xerr = (bin_ends[1:] - bin_ends[:-1]) / 2

    # Now we need to bin our events in order to get an median opening angle (5 degrees/bin)
    binnum = np.digitize(x, bin_ends, right=False) - 1  # Get bin number for each event
    dzen_binned = [[] for i in range(len(bins))]  # create empty array with 1 list for each bin
    for i in range(len(y)):
        dzen_binned[binnum[i]].append(y[i])  # Fill bins with opening angles

    yave = []
    yler = []
    yher = []
    for b in dzen_binned:
        yave.append(np.median(b))
        yler.append(np.quantile(b, .16))
        yher.append(np.quantile(b, .84))

    yerr = [yler, yher]
    return [yave, xerr, yerr, bins]


def zenith_diff_p_fe_by_tier(mc_zen_p, reco_zen_p, mc_zen_fe, reco_zen_fe,
                             bin_ends, titles, cuts):
    data = []
    for (reco_prot_zen, true_prot_zen, reco_fe_zen, true_fe_zen, bin_end, cut) \
            in zip(reco_zen_p, mc_zen_p, reco_zen_fe, mc_zen_fe, bin_ends, cuts):
        proton = zenith_diff(true_prot_zen, reco_prot_zen, bin_end, cut)
        proton.append("Proton")
        proton.append("red")

        iron = zenith_diff(true_fe_zen, reco_fe_zen, bin_end, cut)
        iron.append("Iron")
        iron.append("blue")
        data.append([proton, iron])

    plot_vs_zenith_bytier(data, titles, "Zenith Difference")


def open_angle_vs_zen(true_zen, reco_zen, true_az, reco_az, bin_ends, cut=65):
    opening_angle = find_opening_angle(true_zen, reco_zen, true_az, reco_az)
    # The simulation does not include any zenith values above 65 degrees, so we will cut out values that get too high
    y = np.delete(opening_angle, np.argwhere(reco_zen >= cut))
    x = np.delete(reco_zen, np.argwhere(reco_zen >= cut))

    # Get center of bins
    bins = (bin_ends[1:] + bin_ends[:-1]) / 2
    xerr = (bin_ends[1:] - bin_ends[:-1]) / 2

    # Now we need to bin our events in order to get an median opening angle (5 degrees/bin)
    binnum = np.digitize(x, bin_ends, right=False) - 1  # Get bin number for each event
    dzen_binned = [[] for i in range(len(bins))]  # create empty array with 1 list for each bin
    for i in range(len(y)):
        dzen_binned[binnum[i]].append(y[i])  # Fill bins with opening angles

    yave = []
    yler = []
    yher = []
    for b in dzen_binned:
        yave.append(np.median(b))
        yler.append(np.quantile(b, .16))
        yher.append(np.quantile(b, .84))

    yerr = [yler, yher]
    return [yave, xerr, yerr, bins]


def open_angle_p_fe_by_tier(mc_zen_p, reco_zen_p, mc_az_p, reco_az_p, mc_zen_fe,
                            reco_zen_fe, mc_az_fe, reco_az_fe, bin_ends, titles,
                            cuts):
    data = []

    for (reco_prot_zen, true_prot_zen, reco_prot_az, true_prot_az, reco_fe_zen,
         true_fe_zen, reco_fe_az, true_fe_az, bin_end, cut) in zip(reco_zen_p,
                                                                   mc_zen_p,
                                                              reco_az_p, mc_az_p,
                                                              reco_zen_fe, mc_zen_fe,
                                                              reco_az_fe, mc_az_fe,
                                                              bin_ends, cuts):
        proton = open_angle_vs_zen(true_prot_zen, reco_prot_zen, true_prot_az,
                                   reco_prot_az, bin_end, cut)
        proton.append("Proton")
        proton.append("red")

        iron = open_angle_vs_zen(true_fe_zen, reco_fe_zen, true_fe_az,
                                 reco_fe_az, bin_end, cut)
        iron.append("Iron")
        iron.append("blue")
        data.append([proton, iron])

    plot_vs_zenith_bytier(data, titles, "Opening Angle")

def open_angle_p_fe_first_slice(mc_zen_p, reco_zen_p, mc_az_p, reco_az_p, mc_zen_fe, reco_zen_fe, mc_az_fe, reco_az_fe, binends, reco_name, cuts):
    """ This is more hard coded (2 by 2 as opposed to length of data list. Should be updated to work for variable number of tiers."""
    fig = plt.figure(figsize=(16, 8), constrained_layout=True)
    ax_array = fig.subplots(2, 2, squeeze=False)
    fig.suptitle(r'Opening Angle Distribution for 0 to $16\degree$ Zenith for Each Tier', fontsize = 15)

    x_ind = 0
    y_ind = 0

    ranges = [15, 6, 3, 3]
    for (reco_prot_zen, true_prot_zen, reco_prot_az, true_prot_az, reco_fe_zen, 
         true_fe_zen, reco_fe_az, true_fe_az, bin_ends, reconstruction, ran, cut) in zip(reco_zen_p, mc_zen_p, reco_az_p,
                                                                               mc_az_p, reco_zen_fe, mc_zen_fe, 
                                                                            reco_az_fe, mc_az_fe, binends, reco_name, ranges, cuts):
        dlamda_p = np.deg2rad(reco_prot_az-true_prot_az)
        phi1_p = np.deg2rad(90-reco_prot_zen)
        phi2_p = np.deg2rad(90-true_prot_zen)
        opening_angle_p = np.rad2deg(np.arccos(np.sin(phi1_p)*np.sin(phi2_p) + np.cos(phi1_p)*np.cos(phi2_p)*np.cos(dlamda_p)))
    
        # The simulation does not include any zenith values above 65 degrees, so we will cut out values that get too high
        y_p = np.delete(opening_angle_p, np.argwhere(reco_prot_zen>=cut))
        x_p = np.delete(reco_prot_zen, np.argwhere(reco_prot_zen>=cut))
    
        #bins = bins
        #Get center of bins
        bins  = (bin_ends[1:] + bin_ends[:-1])/2
        xerr = (bin_ends[1:] - bin_ends[:-1])/2
    
        # Now we need to bin our events in order to get an median opening angle (5 degrees/bin)
        binnum_p = np.digitize(x_p, bin_ends, right=False) - 1 # Get bin number for each event 
        dzen_binned_p = [[] for i in range(len(bins))] # create empty array with 1 list for each bin 
        for i in range(len(y_p)):
            dzen_binned_p[binnum_p[i]].append(y_p[i]) # Fill bins with opening angles
        
        yave_p = []
        yler_p = []
        yher_p = []
        for b in dzen_binned_p:
            yave_p.append(np.median(b))
            yler_p.append(np.quantile(b, .16))
            yher_p.append(np.quantile(b, .84))
        
        yerr_p=[yler_p[0], yher_p[0]]
        proton_slice1 = dzen_binned_p[0]
        
        dlamda_fe = np.deg2rad(reco_fe_az-true_fe_az)
        phi1_fe = np.deg2rad(90-reco_fe_zen)
        phi2_fe = np.deg2rad(90-true_fe_zen)
        opening_angle_fe = np.rad2deg(np.arccos(np.sin(phi1_fe)*np.sin(phi2_fe) + np.cos(phi1_fe)*np.cos(phi2_fe)*np.cos(dlamda_fe)))
    
        # The simulation does not include any zenith values above 65 degrees, so we will cut out values that get too high
        y_fe = np.delete(opening_angle_fe, np.argwhere(reco_fe_zen>=cut))
        x_fe = np.delete(reco_fe_zen, np.argwhere(reco_fe_zen>=cut))
        #bins = bins
        #Get center of bins
        #bins  = (bin_ends[1:] + bin_ends[:-1])/2
        #xerr = (bin_ends[1:] - bin_ends[:-1])/2
    
        # Now we need to bin our events in order to get an median opening angle (5 degrees/bin)
        binnum = np.digitize(x_fe, bin_ends, right=False) - 1 # Get bin number for each event 
        dzen_binned_fe = [[] for i in range(len(bins))] # create empty array with 1 list for each bin 
        for i in range(len(y_fe)):
            dzen_binned_fe[binnum[i]].append(y_fe[i]) # Fill bins with opening angles
        
        yave_fe = []
        yler_fe = []
        yher_fe = []
        for b in dzen_binned_fe:
            yave_fe.append(np.median(b))
            yler_fe.append(np.quantile(b, .16))
            yher_fe.append(np.quantile(b, .84))

        yerr_fe=[yler_fe[0], yher_fe[0]]
        iron_slice1 = dzen_binned_fe[0]
        ax_array[x_ind, y_ind].hist(proton_slice1, bins=int(ran/15*100), range=[0,ran], label='Proton', color='red', histtype='step')
        ax_array[x_ind, y_ind].hist(iron_slice1, bins=int(ran/15*100), range=[0,ran], label='Iron', color='blue', histtype='step')
        
        ax_array[x_ind, y_ind].vlines(yave_p[0], 0, 1, transform=ax_array[x_ind, y_ind].get_xaxis_transform(), linestyle='solid', colors='r')
        ax_array[x_ind, y_ind].vlines(yave_fe[0], 0, 1, transform=ax_array[x_ind, y_ind].get_xaxis_transform(), linestyle='solid', colors='b')
        
        ax_array[x_ind, y_ind].vlines(yerr_p, 0, 1, transform=ax_array[x_ind, y_ind].get_xaxis_transform(), linestyle='dashed', colors='r')
        ax_array[x_ind, y_ind].vlines(yerr_fe, 0, 1, transform=ax_array[x_ind, y_ind].get_xaxis_transform(), linestyle='dashed', colors='b')
        
        ax_array[x_ind, y_ind].legend(loc='lower right')
        ax_array[x_ind, y_ind].set_ylabel('Count')
        ax_array[x_ind, y_ind].set_xlabel(r'$Opening\;angle\;[\degree]$ ' + reconstruction)
        #ax_array[x_ind, y_ind].set_ylim()
        ax_array[x_ind, y_ind].set_xlim(0, ran)
        ax_array[x_ind, y_ind].set_xticks(np.arange(0, ran, ran/10))
        locs = ax_array[x_ind, y_ind].get_yticks()
        ax_array[x_ind, y_ind].set_yticklabels(np.round(locs/len(proton_slice1), 2))
        #ax_array[x_ind, y_ind].set_yscale('log')
                
        y_ind = y_ind + 1
        if(y_ind == 2):
            y_ind = 0
            x_ind = x_ind + 1


"""
Weighted quantile
From Stack Overflow: @Alleo
https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy

Function used to find median & error bars of weighted distribuiton 
"""


def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """
    Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!

    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.

    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

def bin_weighted_param(param, bin_ends, param_weights=None):
    # Define lists
    indices = []
    binned = []
    binned_w = []
    err = []
    """
    np.digitize(a, bins) gives the bin number (index+1) that each entry in a should go into
    So by subtracting 1, we can get the index, and create a 2D array with len(zenbins) rows and 
    N_i events in each row (columns)

    e.g. mczenbins[0] will have all the MCPrimary zenith readings between 0-3 degrees
    and brnazbins[3] will have all the BurnSample azimuth readings between 9-12 degrees
    """
    # First, get the indexies of where all zenith events should go using digitize
    indices = np.digitize(param, bin_ends, right=False) - 1
    for b in range(len(bin_ends) - 1):
        binned.append(param[indices == b])
        if param_weights is not None:
            binned_w.append(param_weights[indices == b])

    if len(binned_w) > 0:
        for b in binned_w:
            err.append(math.sqrt(np.sum(b * b)))
    else:
        for b in binned:
            err.append(math.sqrt(len(b)))
    return [binned, binned_w, err]


def bin_particle_energy(particle_param, particle_en, particle_weights, bin_ends):
    median = []
    err = []
    indices = np.digitize(particle_param, bin_ends, right=False) - 1
    for b in range(len(bin_ends) - 1):
        m, u, l = weighted_quantile(particle_en[indices == b] * (10 ** -6), [.5, .84, .16],
                                    sample_weight=particle_weights[indices == b])
        median.append(m)
        err.append([m - l, u - m])
    return [median, err]

def compare_dist(reco, reco_w, data, reco_name, param_name, bin_ends):
    reco_bin, reco_bin_w, reco_err = bin_weighted_param(reco, bin_ends, reco_w)
    data_bin, data_bin_w, data_err = bin_weighted_param(data, bin_ends)

    bin_centres = (bin_ends[1:] + bin_ends[:-1])/2

    y_reco = plt.hist(reco, bin_ends, weights=reco_w)[0]
    y_data = plt.hist(data, bin_ends)[0]
    plt.close()

    fig = plt.figure(figsize=(10, 7), constrained_layout=True)
    ax_array = fig.subplots(2, 2, squeeze=False)
    fig.suptitle('Sim vs Data '+ reco_name + ' ' + param_name)

    ax_array[0,0].set_title('Linear')
    ax_array[0,0].set_xlabel(r'$Zenith\;angle\;[\degree]$')
    ax_array[0,0].set_ylabel('Fraction of Events')
    ax_array[0,0].errorbar(bin_centres, y_reco/sum(y_reco), yerr=reco_err/sum(y_reco), label='Sim ' + reco_name, capsize=3)
    ax_array[0,0].errorbar(bin_centres, y_data/sum(y_data), yerr=data_err/sum(y_data), label='Data ' + reco_name, capsize=3)
    ax_array[0,0].set_ylim(0, .1)
    ax_array[0,0].legend(loc='lower center')

    ax_array[0,1].set_title('Log')
    ax_array[0,1].set_xlabel(r'$Zenith\;angle\;[\degree]$')
    ax_array[0,1].set_ylabel('Fraction of Events')
    ax_array[0,1].semilogy()
    ax_array[0,1].errorbar(bin_centres, y_reco/sum(y_reco), yerr=reco_err/sum(y_reco), label='Sim ' + reco_name, capsize=3)
    ax_array[0,1].errorbar(bin_centres, y_data/sum(y_data), yerr=data_err/sum(y_data), label='Data ' + reco_name, capsize=3)
    ax_array[0,1].legend(loc='lower center')


    ax_array[1,0].set_title('Ratio Linear')
    ax_array[1,0].set_xlabel(r'$Zenith\;angle\;[\degree]$')
    ax_array[1,0].set_ylabel('Ratio')
    ax_array[1,0].errorbar(bin_centres, ((y_reco/sum(y_reco))/(y_data/sum(y_data))), yerr=.1, label='Ratio', capsize=3)
    ax_array[1,0].legend(loc='lower center')


    ax_array[1,1].set_title('Ratio Log')
    ax_array[1,1].set_xlabel(r'$Zenith\;angle\;[\degree]$')
    ax_array[1,1].set_ylabel('Ratio')
    ax_array[1,1].set_yscale('symlog')
    ax_array[1,1].errorbar(bin_centres, ((y_reco/sum(y_reco))/(y_data/sum(y_data))), yerr=.1, label='Ratio', capsize=3)
    ax_array[1,1].legend(loc='lower center')


def dist_sim_data_all_tiers(data_list, reco_list, param, bin_ends, titles, logscale=False, output='./plot.png'):
    n = len(data_list)
    height = math.ceil(n / 2)
    fig = plt.figure(figsize=(12, 6 * height), constrained_layout=True)
    ax_array = fig.subplots(2, height, squeeze=False)
    if logscale:
        fig.suptitle('Log Recon. Sim vs Data for ' + param + ' By Tier',
                 fontsize=30)
    else:
        fig.suptitle('Linear Recon. Sim vs Data for ' + param + 'By Tier',
                     fontsize=30)
        
    bin_centres = (bin_ends[1:] + bin_ends[:-1])/2
    
    row = 0
    col = 0

    for data, reco, title in zip(data_list, reco_list,
                                       titles):
        y_reco = np.histogram(reco[0], bin_ends, weights=reco[1])[0]
        y_data = np.histogram(data[0], bin_ends)[0]

        if logscale:
            ax_array[row, col].semilogy()

        ax_array[row, col].set_ylabel('Fraction of Events', fontsize=20)
        ax_array[row, col].errorbar(bin_centres, y_reco/sum(y_reco), yerr=reco[2]/sum(y_reco), label='Sim ' + title[1], capsize=3)
        ax_array[row, col].errorbar(bin_centres, y_data/sum(y_data), yerr=data[1]/sum(y_data), label='Data ' + title[1], capsize=3)
        #ax_array[row, col].set_ylim(0, .1)
    
        ax_array[row, col].set_title('Tier {}'.format(title[0]), fontsize=20)
        ax_array[row, col].set_xlabel(param + ' angle' + '$[\degree]$', fontsize=20)
        ax_array[row, col].legend(loc='lower center', fontsize=15)
    
        col = col + 1
        if col == 2:
            col = 0
            row = row + 1
        
    plt.savefig(output, transparent=True)
            
def ratio_sim_data_all_tiers(data_list, reco_list, param, bin_ends, titles, logscale=False, output='./plot.png'):
    n = len(data_list)
    height = math.ceil(n / 2)
    fig = plt.figure(figsize=(12, 6 * height), constrained_layout=True)
    ax_array = fig.subplots(2, height, squeeze=False)
    if logscale:
        fig.suptitle('Ratio Log Recon. Sim over Data for ' + param + ' By Tier',
                 fontsize=30)
    else:
        fig.suptitle('Ratio Linear Recon. Sim over Data for ' + param + 'By Tier',
                     fontsize=30)
        
    bin_centres = (bin_ends[1:] + bin_ends[:-1])/2
    
    row = 0
    col = 0

    for data, reco, title in zip(data_list, reco_list,
                                       titles):
        y_reco = np.histogram(reco[0], bin_ends, weights=reco[1])[0]
        y_data = np.histogram(data[0], bin_ends)[0]

        if logscale:
            ax_array[row, col].set_yscale('symlog')
            #ax_array[row, col].semilogy()
            ax_array[row, col].set_ylim(0, 10)
            ax_array[row, col].set_yticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
            ax_array[row, col].set_yticklabels(["1", "", "", "", "", "", "", "", "", "10"])
        else:
            ax_array[row, col].set_ylim(0, 2)
        ax_array[row, col].set_title('Tier {}'.format(title[0]), fontsize=20)
        ax_array[row, col].set_xlabel(param + ' angle' + '$[\degree]$', fontsize=20)
        ax_array[row, col].set_ylabel('Ratio', fontsize=20)
        ax_array[row, col].errorbar(bin_centres, ((y_reco/sum(y_reco))/(y_data/sum(y_data))), yerr=.1, label='Ratio', capsize=3)
        ax_array[row, col].legend(loc='lower center', fontsize=20)
    
        col = col + 1
        if col == 2:
            col = 0
            row = row + 1
    plt.savefig(output, transparent=True)
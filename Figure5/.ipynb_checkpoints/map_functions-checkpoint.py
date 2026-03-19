#!/usr/bin/env python

import healpy as hp
import numpy as np

from matplotlib.colors import LinearSegmentedColormap

def adjust_colorbar_limits(skymap, cb_ticks):

    # Take no action if none of the colorbar parameters are set
    if cb_ticks==None:
        return skymap
    
    # Identify masked pixels
    mask = np.ma.getmaskarray(skymap)
    
    # Identify unmasked pixels below colorbar minimum
    under_cut = np.logical_not(mask) * (skymap < cb_ticks[0])
    if under_cut.sum() > 0:
        print(f'Warning: {under_cut.sum()} pixels found with values under your colorbar minimum!')
        print('Setting to colorbar minimum (view them by removing cb_min/cb_max/cb_ticks)')
        skymap[under_cut] = cb_ticks[0]

    # Identify unmasked pixels above colorbar maximum
    above_cut = np.logical_not(mask) * (skymap > cb_ticks[-1])
    if above_cut.sum() > 0:
        print(f'Warning: {above_cut.sum()} pixels found with values over your colorbar maximum!')
        print('Setting to colorbar maximum (view them by removing cb_min/cb_max/cb_ticks)')
        skymap[above_cut] = cb_ticks[-1]

    return skymap
    

def maskMap(skymap, decmin=-90., decmax=-35.):
    
    npix  = len(skymap)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))

    thetaMin = np.radians(90 - decmin)
    thetaMax = np.radians(90 - decmax)
    mask = (theta <= thetaMin) * (theta >= thetaMax)
    
    # Apply the mask
    masked_map = hp.ma(skymap)
    masked_map.mask = np.logical_not(mask).astype(bool)

    return masked_map


def smoothMap(m, smooth_deg, norm=False):

    npix  = len(m)
    nside = hp.npix2nside(npix)
    smooth_map = np.zeros(npix)

    vec = np.transpose(hp.pix2vec(nside, np.arange(npix)))
    for i in range(npix):
        neighbors = hp.query_disc(nside, vec[i], np.radians(smooth_deg))
        smooth_map[i] += m[neighbors].sum()
        if norm:
            smooth_map[i] /= (len(neighbors) + 1)

    return smooth_map


def SetupAbsThresholdColormap(skymap, threshold, cb_ticks=None):
    """ Create a color map for "two-sided" thresholds.  Below the threshold,
        the map is a cool green-blue palette.  Between the lower and upper
        threshold, the map is gray-white-gray.  Above the upper threshold,
        the map is a warm red-yellow palette.
    """

    # Extract minimum and maximum colorbar values from plot if not provided
    pixel_cut = (skymap!=hp.UNSEEN) * (skymap!=np.inf) * (skymap!=-np.inf) * (skymap==skymap)
    unmasked_pix = skymap[pixel_cut]  # Eliminates masked pixels, infinities, and nans
    if cb_ticks == None:
        cb_ticks = [unmasked_pix.min(), unmasked_pix.max()]

    # In the case of a threshold outside of our range, set range to threshold
    if cb_ticks[0] > -threshold:
        cb_ticks[0] = -threshold
    if cb_ticks[-1] < threshold:
        cb_ticks[-1] = threshold

    x1 = (-threshold - cb_ticks[0]) / (cb_ticks[-1] - cb_ticks[0])
    x3 = (cb_ticks[-1] - threshold) / (cb_ticks[-1] - cb_ticks[0])
    x2 = 1. - x1 - x3
    gvl = 0.5
    thrDict = {
        "red"    : ((0.0, 1.0, 0.5), (x1, 0.0, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                    (x1 + x2, gvl, 0.7), (1.0, 1.0, 1.0)),
        "green"  : ((0.0, 1.0, 1.0), (x1, 0.0, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                    (x1 + x2, gvl, 0.0), (1.0, 1.0, 1.0)),
        "blue"   : ((0.0, 1.0, 1.0), (x1, 0.7, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                    (x1 + x2, gvl, 0.0), (1.0, 0.5, 1.0)) }
    
    return LinearSegmentedColormap("thresholdColormap", thrDict, 256)
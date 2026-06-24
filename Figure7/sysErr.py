#!/usr/bin/env python
# coding: utf-8

#imports needed packages
import healpy as hp
import numpy as np
import argparse


if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Makes angular power spectrum error bars')
    p.add_argument('-f', '--infile', dest='infile',
            help='Input filename')
    p.add_argument('-n', '--n', dest='n', 
            type=int,
            default = int(1e3),
            help='Number of simulated power spectra to generate')
    p.add_argument('-o', '--out', dest='out', 
            help='Output filename')

    args = p.parse_args()

    # Load combined relative intensity map
    data, bg, relint = hp.read_map(args.infile, range(3))

    # Map parameters
    npix = relint.size
    nside = hp.npix2nside(npix)
    lmax = 3*nside - 1

    # Masking
    theta, phi = hp.pix2ang(nside, range(npix))
    mask = (theta < np.pi - np.radians(55.))
    relint[mask] = 0
    bg[mask] = 0

    # Background map used for weighting
    weight = bg / bg[~mask].mean()

    # Apply weight and shift to keep average relative intensity at 0
    relint_w = relint * weight
    relint_w -= np.average(relint_w)

    # Calculate power spectrum
    Cl = hp.anafast(relint_w, lmax=lmax)

    # Generate fake maps, calculate power spectra, and store
    fakeCl = np.zeros((args.n+1, lmax+1))
    for n in range(args.n):
        fakeRL = hp.synfast(Cl, nside=nside)
        # Calculate and store power spectrum
        fakeCl[n] = hp.anafast(fakeRL, lmax=lmax)

    # Add real Cl at the end of fakeCl for storage
    fakeCl[-1] = Cl


    # Save fakeCl
    np.save(f'{args.out}.npy', fakeCl)


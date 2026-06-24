#!/usr/bin/env python
# coding: utf-8

import healpy as hp
import numpy as np
import argparse


if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Makes angular power spectrum isotropic bands')
    p.add_argument('-f', '--infile', dest='infile',
            help='Input filename')
    p.add_argument('-n', '--n', dest='n',
            type=int, default=int(1e3),
            help='Number of simulated power spectra to generate')
    p.add_argument('-o', '--out', dest='out',
            help='Output filename')

    args = p.parse_args()

    # Background map should represent detector response to isotropic sky
    data, bg, _ = hp.read_map(args.infile, range(3))

    # Map paramaters
    npix = bg.size
    nside = hp.npix2nside(npix)
    lmax = 3*nside - 1

    # Masking
    theta, phi = hp.pix2ang(nside, range(npix))
    mask = (theta < np.pi - np.radians(55.))
    bg[mask] = 0

    # "Normalize" by dividing by average value of unmasked pixels
    weight = bg / bg[~mask].mean()

    # Generate fake maps, calculate power spectra, and store
    fakeCl = np.zeros((args.n, lmax+1))

    for n in range(args.n):

        # Wiggle isotropic map within poisson uncertainties
        dummyMap = np.random.poisson(bg)

        # Calculate relative intensity
        relint = (dummyMap/bg - 1)
        relint[relint != relint] = 0    # deal with NaN's

        # Weight, shift, calculate power spectrum
        relint_w = relint * weight
        relint_w -= np.average(relint_w)
        fakeCl[n] = hp.anafast(relint_w, lmax=lmax)


    # Save fakeCl
    np.save(f'{args.out}.npy', fakeCl)
    print(f'Isotropic Cl values written to {args.out}.npy')



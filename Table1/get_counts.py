#!/usr/bin/env python

import numpy as np
import healpy as hp
import argparse

from pathlib import Path
from datetime import datetime as dt


def get_counts(tier, year, outfile):

    # Load in maps
    root = Path('/data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output')
    root = root / 'sidereal_unblinded' / f'tier{tier}' / 'fitsbydate'

    # Collect all files associated with a given tier
    print(f'Collecting map files for Tier {tier}...')
    files = np.array(sorted(root.glob('????-??-??/*.fits.gz')))

    # Also find the date for each
    dates = np.array([dt.strptime(f.parent.name, '%Y-%m-%d') for f in files])
    print(f'  {len(files)} files found')

    # Establish IceCube calendar years
    start = dt(year, 5, 13)
    stop  = dt(year+1, 5, 13)
    cut   = (dates >= start) * (dates < stop)
    print(f'\nReading in {cut.sum()} files for {year}...')

    # Extract count information from relevant files
    count = 0
    for year_file in files[cut]:
        data = hp.read_map(year_file)
        count += data.sum()

    print(f'\nWriting results to {outfile}...\n')
    with open(outfile, 'w') as f:
        f.write(str(count))


if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Calculates event counts for IceTop 11-year paper')
    p.add_argument('-t', '--tier', 
            dest='tier', type=int,
            choices=[1,2,3,4],
            help='Energy tier selection')
    p.add_argument('-y', '--year', 
            dest='year', type=int,
            choices=[yy for yy in range(2011,2022)],
            help='IceTop calendar year (begins on May 13)')
    p.add_argument('-o', '--outfile', 
            dest='outfile', type=str,
            help='Destination for output')
    args = p.parse_args()

    get_counts(args.tier, args.year, args.outfile)



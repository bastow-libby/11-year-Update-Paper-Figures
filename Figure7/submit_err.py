#!/usr/bin/env python

import argparse
from pathlib import Path
from submitter.pysubmit import pysubmit

if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Submits jobs for power spec uncertainty calculations')
    p.add_argument('-n', '--n', dest='n',
            type=int, default=int(1e3),
            help='Number of simulated power spectra to generate per job')
    p.add_argument('-N', '--N', dest='N',
            type=int, default=int(100),
            help='Number of jobs to submit')
    p.add_argument('--err', dest='err',
            choices=['iso','stat','sys'],
            help='Select which uncertainties to calculate [iso|stat|sys]')
    p.add_argument('-o', '--out', dest='out',
            default='/data/user/fmcnally/it_anisotropy/powerspec',
            help='Output directory')

    args = p.parse_args()


    tiers = [i for i in range(1, 5)]

    # Location for data storage
    out_path = Path(args.out)
    out_path.mkdir(parents=True, exist_ok=True)

    # Directory for this file
    cwd = Path(__file__).resolve().parent

    # Header establishes virtual environment for submission
    cvmfs = '/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/setup.sh'
    pyenv = '/home/fmcnally/i3/442/bin/activate'
    header = ['#!/bin/bash', f'eval $({cvmfs})', f'source {pyenv}']

    for tier in tiers:

        # Final combined relative intensity map
        map_dir  = Path('/data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output')
        tier_dir = map_dir / f'sidereal_unblinded/tier{tier}/reconstruction/'
        map_files = tier_dir.glob('relintensityiter/combined*.fits.gz')
        map_file = sorted(map_files)[-1]

        for i in range(args.N):

            out = f'{args.out}/t{tier}_{args.err}_{i:03}'
            cmd = f'{cwd}/{args.err}Err.py -f {map_file} -n {args.n} -o {out}'

            # Submit to cluster
            jobID = f'{Path(out).name}'
            pysubmit(cmd, jobID=jobID, header=header)



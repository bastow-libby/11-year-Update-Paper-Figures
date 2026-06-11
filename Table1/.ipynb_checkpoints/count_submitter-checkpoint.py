#!/usr/bin/env python

from pathlib import Path
from submitter.pysubmit import pysubmit

if __name__ == "__main__":

    tiers = [i for i in range(1, 5)]
    years = [yy for yy in range(2011, 2022)]

    # Location for data storage
    out_base = Path('/data/user/fmcnally/it_anisotropy/event_counts')
    out_base.mkdir(parents=True, exist_ok=True)

    for tier in tiers:

        for year in years:

            # Check to see if output file already exists
            out_file = out_base / f'tier{tier}_{year}.txt'
            if out_file.is_file():
                print(f'File {out_file.name} already exists. Skipping...')
                continue

            # Executable with arguments
            cwd = Path(__file__).parent.resolve()
            ex = f'{cwd}/get_counts.py -t {tier} -y {year} -o {out_file}'

            # Header establishes virtual environment
            cvmfs = '/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/setup.sh'
            pyenv = '/home/fmcnally/i3/442/bin/activate'
            header = ['#!/bin/bash', f'eval $({cvmfs})', f'source {pyenv}']

            # Memory request required for higher tiers (not sure why...)
            sublines = ['request_memory = 4000']

            # Submit to cluster
            jobID = f'counts_{out_file.name}'
            pysubmit(ex, jobID=jobID, header=header, sublines=sublines)

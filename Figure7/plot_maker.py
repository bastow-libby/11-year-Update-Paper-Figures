#!/usr/bin/env python

import subprocess
from pathlib import Path

def plot_aps(map_file, out_dir, err_dir, tier, label, 
        iso_label=False, icecube_prelim=False):

    # Code to make Angular Power Spectrum
    cwd = Path.cwd()
    cmd = f'{cwd}/aps.py -f {map_file}'

    # Add uncertainty files to the plot (if present)
    iso_file = Path(f'{err_dir}/t{tier}_iso.npy')
    if iso_file.is_file():
        cmd += f' -i {iso_file}'
    else: print('Isotropic error bands not found!')

    stat_file = Path(f'{err_dir}/t{tier}_stat.txt')
    if stat_file.is_file():
        cmd += f' -st {stat_file}'
    else: print('Statistical uncertainties not found!')

    sys_file = Path(f'{err_dir}/t{tier}_sys.txt')
    if sys_file.is_file():
        cmd += f' -sy {sys_file}'
    else: print('Systematic uncertainties not found!')

    # Labeling
    if iso_label:
        cmd += f' -il'
    if icecube_prelim:
        cmd += f' -icp'

    # Output name for image
    cmd += f' -o {out_dir}/APS_T{tier}.pdf -l {label}'
    subprocess.run(cmd.split(' '))


if __name__ == "__main__":

    '''
    NOTE: the code originally had the following pre-compiled metaproject
    as a header. It's been removed as it doesn't have support for healpy
    #!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/icetray-start
    #METAPROJECT: icetray/v1.17.0
    '''

    # File paths
    map_dir  = Path('/data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output')
    out_dir = '.'
    err_dir = '/data/user/fmcnally/it_anisotropy/powerspec'

    # Energy labels
    labels = {1:'280_TeV', 2:'900_TeV', 3:'2300_TeV', 4:'6600_TeV'}

    for tier, label in labels.items():

        tier_dir = map_dir / f'sidereal_unblinded/tier{tier}/reconstruction/'
        map_files = tier_dir.glob('relintensityiter/combined*.fits.gz')
        map_file = sorted(map_files)[-1]

        # Suppress iso label on tiers 2-4
        iso_label = True if tier==1 else False

        plot_aps(map_file, out_dir, err_dir, tier, label, iso_label=iso_label)

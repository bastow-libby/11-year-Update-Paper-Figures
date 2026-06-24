#!/usr/bin/env python

import numpy as np
from pathlib import Path
import argparse


def merge_iso_files(file_list, out):

    # Load numpy array files
    fake_Cl = []
    for iso_file in file_list:
        fCl_i = np.load(iso_file)
        fake_Cl.append(fCl_i)

    # Collapse to a single array
    fake_Cl = np.concatenate(fake_Cl, axis=0)

    # Calculate 1-, 2-, and 3-sigma containment values from simulated Cls
    lims = [0.15, 2.5, 16, 84, 97.5, 99.85]
    fCl = np.percentile(fake_Cl, lims, axis=0)

    # Remove monopole (l=0)
    fCl = fCl[:,1:]

    iso = {}
    vals = [-3, -2, -1, 1, 2, 3]
    for i, val in enumerate(vals):
        iso[val] = fCl[i]

    # Save as dictionary
    np.save(out, iso)
    print(f'Isotropic noise bands written to {Path(out).name}')


def merge_stat_files(file_list, out):

    # Load numpy array files
    fake_Cl = []
    for stat_file in file_list:
        fCl_i = np.load(stat_file)
        fake_Cl.append(fCl_i)

    # Collapse to a single array
    fake_Cl = np.concatenate(fake_Cl, axis=0)

    # Calculate 1-sigma containment values from simulated Cls
    lims = [16, 50, 84]
    s_16, s_50, s_84 = np.percentile(fake_Cl, lims, axis=0)
    dCl = [ s_50-s_16, s_84-s_50 ]

    # Save as dictionary
    np.savetxt(out, dCl)
    print(f'Statistical uncertainties written to {Path(out).name}')


def merge_sys_files(file_list, out):

    # Load numpy array files
    fake_Cl = []
    for sys_file in file_list:
        fCl_i = np.load(sys_file)
        Cl = fCl_i[-1]      # Last entry reserved for real Cl values
        fake_Cl.append(fCl_i[:-1])

    # Collapse to a single array
    fake_Cl = np.concatenate(fake_Cl, axis=0)

    # Calculate 1-sigma containment values from simulated Cls
    lims = [16, 84]
    fCl = np.percentile(fake_Cl, lims, axis=0)

    # Convert 1-sigma containments to errorbar format for plotting
    dCl = [ Cl-fCl[0], fCl[1]-Cl ]

    # Save to text file
    np.savetxt(out, dCl)
    print(f'Systematic uncertainties written to {Path(out).name}')


if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description=('Calculates isotropic, systematic, and statistical '
                         'uncertainties from isoErr, sysErr, and statErr '
                         'output files'))
    p.add_argument('--file_dir', dest='file_dir',
            #default=('/data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output'
            #         '/icecube/powerspec')
            default='/data/user/fmcnally/it_anisotropy/powerspec',
            help='Directory with output files')
    p.add_argument('-o', '--out', dest='out',
            #default=('/data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output'
            #         '/icecube/powerspec')
            default='/data/user/fmcnally/it_anisotropy/powerspec',
            help='Output directory (defaults to same as input)')

    args = p.parse_args()

    print(f'Output directory: {args.out}')

    for tier in range(1, 5):

        print(f'Merging files for tier {tier}...')

        iso_files = Path(args.file_dir).glob(f't{tier}_iso_???.npy')
        out = f'{args.out}/t{tier}_iso.npy'
        merge_iso_files(iso_files, out)

        stat_files = Path(args.file_dir).glob(f't{tier}_stat_???.npy')
        out = f'{args.out}/t{tier}_stat.txt'
        merge_stat_files(stat_files, out)

        sys_files = Path(args.file_dir).glob(f't{tier}_sys_???.npy')
        out = f'{args.out}/t{tier}_sys.txt'
        merge_sys_files(sys_files, out)


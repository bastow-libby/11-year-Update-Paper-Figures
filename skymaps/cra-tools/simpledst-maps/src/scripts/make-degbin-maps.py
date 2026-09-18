#!/bin/env python

import re, argparse
import random
from pathlib import Path
import getpass

from npx4.pysubmit import pysubmit


if __name__ == "__main__":

    user = getpass.getuser()

    p = argparse.ArgumentParser(
            description=('Makes degree-binned healpix maps for use with '
                         'Ahlers likelihood background estimation method'))

    # Common options
    p.add_argument('-y', '--year', dest='year',
            nargs='+',
            help='Detector year [2011-2021]')
    p.add_argument('--outdir', dest='outdir',
        default='/data/user/@USER_DIR@/burnsample',
            help='Destination directory for output files')
    p.add_argument('--outfile', dest='outfile',
            help='Base name for output file')
    p.add_argument('-m', '--method', dest='method',
            choices=['sid','solar','anti','ext'],
            default='sid',
            help='Time frame for coordinate conversion')

    # Additional options
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Option for running off cluster to test')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Option to overwrite existing map files')

    # Uncommon options -- use with caution
    p.add_argument('--N', dest='N',
            default=50, # Job length of ~1 hour
            help='Batch size for job submission (# days/job)')
    p.add_argument('--smin', dest='smin',
            default=3.0,
            help='Minimum N_stations, inclusive')
    p.add_argument('--smax', dest='smax',
            default=100.0,
            help='Maximum N_stations, non-inclusive')
    p.add_argument('--dirreco', dest='dirreco',
            choices=['ShowerPlane','Laputop','LaputopSmall'],
            default='ShowerPlane',
            help='Directional reconstruction to use')
    p.add_argument('--sd', dest='sundp',
            default=False, action='store_true',
            help='Correct for solar dipole on event-by-event basis')
    p.add_argument('--nside', dest='nsideout',
            default=64,
            help='Healpix NSide parameter')

    # Temporarily deprecated options
    #p.add_argument('--submit_dir', dest='submit_dir',
    #        default='@CRA_BUILD@/submit-dir',
    #        help='submit directory')
    #p.add_argument('--cvmfs', dest='cvmfs',
    #        default='@CVMFS_SROOTBASE@',
    #        help='CVMFS environment')
    #p.add_argument('--condor-priority', dest='priority'
    #        default=1, 
    #        help='Job priority in condor')

    args = p.parse_args()

    # Default base naming scheme for outfiles (if not provided)
    if args.outfile == None:
        args.outfile = f'ITpass2_{args.method}'
        if args.sundp:
            args.outfile += '_sd'

    # Options to be passed directly to make-local-maps.cc
    c_opts = {'config':'ITpass2'}
    for key in ['smin','smax','method','dirreco','nsideout']:
        c_opts[key] =  vars(args)[key]

    # Print options for user
    print('Parameters for submission:')
    for key, value in sorted(c_opts.items()):
        print('  --%s %s' % (key, value))
    if args.sundp:
        print('  --sundp')

    # Environment for script
    cvmfs = '/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh'
    pyenv = '/home/@USER_DIR@/venv/bin/activate'
    header = ['#!/bin/bash', f'eval $({cvmfs})', f'source {pyenv}']

    # Memory requirements appear to cap out around 2500 MB
    sublines = ["request_memory = 3000"]

    # Lazy approach for all years
    if args.year == ['all']:
        args.year = list(range(2011,2022))
    if args.year == ['tier1']:
        args.year = list(range(2011, 2015))

    # Get working directory for this script
    script_path = Path(__file__).resolve()
    script_dir = script_path.parent

    for year in args.year:

        print(f'Working on {year}...')

        # Collect input files
        fpath = '/data/user/@USER_DIR@/burnsample_root'
        file_list = sorted(Path(fpath).glob(f'*_IC86.{year}_*.root'))

        # Loop over file list to build a list of files for submission
        c_list = []

        for root_file in file_list:

            # Check for existing files
            date = re.findall(r'\d{4}_\d{4}\.root', root_file.name)[-1][:-5]
            yyyy, mmdd = date.split('_')
            test_dir = Path(f'{args.outdir}/IC86.{year}/{yyyy}/{mmdd}')
            test_base = f'CR_ICETOP_LOCAL_3-100S_NSIDE64'

            if not args.overwrite:

                if test_dir.is_dir():

                    # Get list of existing filenames
                    ex_files = {f.name for f in 
                            test_dir.glob(f'{test_base}_degbin-*.fits.gz')}

                    # Generate the string names we expect to see
                    ex_names = (f'{test_base}_degbin-{i:03}.fits.gz' 
                            for i in range(360))

                    # Optimization: fast set comparison
                    if not any(n not in ex_files for n in ex_names):
                        continue
            
            # Create output directory if it doesn't already exist
            test_dir.mkdir(parents=True, exist_ok=True)

            # Remove existing files in case of overwrite or partial coverage
            test_files = [test_dir / f'{test_base}_degbin-{i:03}.fits.gz'
                          for i in range(360)]
            for test_file in test_files:
                test_file.unlink(missing_ok=True)

            # Convert c_opts dictionary to a string for submission
            c_opts['outdir'] = test_dir
            c_opts['outfile'] = test_base
            c_opts['input'] = str(root_file)
            c_str = ' '.join(f'--{k} {v}' for k, v in sorted(c_opts.items()))

            # Include additional flags without arguments at end
            if args.sundp:
                c_str += ' --sundp'

            # Save c_opts for days that pass above checks
            c_list += [c_str]        


        # Batch the c_opts
        sub_lists = [c_list[i:i+args.N] for i in range(0, len(c_list), args.N)]

        # Limit submissions if running in test mode
        if args.test:
            sub_lists = sub_lists[:1]

        # Submit
        cmd = f'/home/@USER_DIR@/11-year-Update-loyola/skymaps/cra-tools/simpledst-maps/build/bin/make-degbin-maps'
        #cmd='/home/fmcnally/cra-tools/simpledst-maps/build/bin/make-local-maps'
        #cmd  = '@CRA_BUILD@/bin/make-local-maps'
        for sub_batch in sub_lists:

            # Name job
            jobID = f'{args.outfile}'
            for s in [sub_batch[0], sub_batch[-1]]:
                date = re.findall(r'IC86\.\d{4}_\d{4}_\d{4}', s)[-1]
                jobID += f'_{date}'

            # Create and print list of executables
            print(f'\n{jobID}')
            ex_list = []
            for c_opts_str in sub_batch:
                ex = f'{cmd} {c_opts_str}'
                ex_list += [ex]
                print(ex)

            pysubmit(ex_list, sublines=sublines, test=args.test, jobID=jobID,
                    header=header
                    #outdir=submit_dir,
                    #priority=args.priority
                    )


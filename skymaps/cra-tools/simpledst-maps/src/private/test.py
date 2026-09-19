#!/usr/bin/env python

import subprocess
from pathlib import Path

if __name__ == "__main__":

    cmd = '/home/fmcnally/cra-tools/simpledst-maps/build/bin/TimeScramble'
    fpath = '/data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/unblinded_root'
    file_list = sorted(Path(fpath).glob('*_IC86.2011_*.root'))
    file_list = [str(f) for f in file_list]
    file_list = file_list[:1]


    ex = f'{cmd} --help'

    file_str = ' '.join(file_list)
    ebins = [4, 4.25, 4.5]
    ebin_str = ' '.join([str(i) for i in ebins])

    spline = '/data/ana/CosmicRay/Anisotropy/IceCube/twelve_year/sim'
    spline = f'{spline}/IC86_20904_hist_spline.fits'

    outdir = '/home/fmcnally/cra-tools/simpledst-maps/src/private'
    outfile = 'TS_test'

    ex = f'{cmd} --input {file_str} --detector IT --outdir {outdir}'
    ex = f'{ex} --outfile {outfile}'
    print(ex)

    subprocess.run(ex.split(' '))

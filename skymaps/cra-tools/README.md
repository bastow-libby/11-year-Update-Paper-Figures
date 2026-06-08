# cra-tools

This is a collection of tools used for analyzing the angular distribution of cosmic-ray arrival directions.


The project is divided into three subprojects:

## TimeScramble
Time-scrambling code for IceCube CRA studies. Currently only reads ROOT format SimpleDST files.

## iter-lhreco
C++ implementation of maximum-likelihood technique for reconstructing cosmic-ray anisotropy maps
http://iopscience.iop.org/article/10.3847/0004-637X/823/1/10

## simpledst-maps
C++ code for extracting data from IceCube's reduced data format (simple-dst) in ROOT/HDF5 and generating local HEALpix maps.


## scripts
Collection of scripts for driving production of extraction, reconstruction and analysis of cosmic-ray data.


## Installation


**Prerequisites**:

iter-lhreco dependencies:

* CFITSIO: depends on cfitsio: https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html 
* GNU Scientific Library: depends on gsl: https://www.gnu.org/software/gsl/doc/html/index.html
* HEALpix: depends on healpix: http://healpix.sourceforge.net/
* BOOST C++ libraries: https://www.boost.org/


simpledst-maps dependencies:

* CFITSIO: depends on cfitsio: https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html 
* GNU Scientific Library: depends on gsl: https://www.gnu.org/software/gsl/doc/html/index.html
* Photospline: extraction depends on photospline to make energy cuts: https://github.com/IceCubeOpenSource/photospline
* Starlink Positional Astronomy Library: extraction depends on pal: https://github.com/Starlink/pal
* CERN ROOT: extraction depends on root: https://root.cern.ch/
* HEALpix: depends on healpix: http://healpix.sourceforge.net/

You will also need:

* [CMake](https://cmake.org) >= 3.1
* A C++11-compliant compiler (e.g. [gcc](https://gcc.gnu.org) >= 4.8.1 or [clang](https://clang.llvm.org) >= 3.3)


**Installation**:

The C++ projects are built using CMake.
Each project has a CMakeList.txt that will detect dependencies and generate a MakeFile. To do this cd into the build directory and excecute the commands:

  cmake ../src;
  make

## Instructions
**Step 1. Build the Project with CMake**:
Run cvmfs (make sure you run the same version you will hardcode into the scripts make-local-maps-IT.py.in and merge-reco.py.in, as of June 2026, that version is 4.3.0)
cd 11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps
mkdir build
cd build
cmake ../src
make

**Step 2. Edit Code for the User**:
Navigate to 11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/scripts. We will edit 2 scripts here.

The first is make-local-maps-IT.py.in. This makes the local maps. Edit these lines:

On line 19, edit the fpath to the user directory where you have stored the burnsample/full dataset root files.
   fpath = '/data/.../root'
On line 65, edit the output directory to where you want the output directory to be (recommend /data/user/@USER_DIR@/... for burnsample, /data/ana/... for the full dataset):
   default='/data/.../OUTPUT\_DIRECTORY'
On line 77, edit the line to the version of cvmfs used in cmake:
    default='/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0',
On line 81, edit the submit directory to the path where your submit directory is:
    default='/home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/submit'
On line 223, update the cmd line to the location of your 'make-local-maps' file:
    cmd='/home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/build/bin/make-local-maps'

The second is merge-reco.py.in. This script combines the local maps and uses the iteration method to estimate relative intensity and significance. Edit these lines:

On line 61, update the line to your output directory:
   default='/data/.../OUTPUT\_DIRECTORY',
On line 73, update the line to the version of cvmfs used in cmake:
   default='/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0',
On line 76, update the line to the path to your submit directory:
   default='/home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/submit',
On line 214, update the line to the location of your 'combine-local-maps' folder:
   cmd  = '/home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/build/bin/combine-local-maps'
On line 226, update the line to:
   cmd  = '/home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/build/bin/illh-reco'
On line 228, update the line to:
   '--cfg','/home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/build/resources/{detector}.json'.format(\*\*s\_opts),

We also have to edit 11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/private/SimpleDST.cc depending on which energy tier we intend to process:

For Tier 1, comment out line 39 and uncomment line 40, so it reads
    reco = "ShowerPlane";
For Tiers 2, 3 or 4, uncomment out line 39, and comment line 40, so it reads
    reco = "Laputop";

**Step 3. Rebuild the Project**
If you have edited a .cc file, remake the project as follows:
run cvmfs (use the same version used in all previous steps)
cd 11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps
cd build
cmake ../src
make

**Step 4. Edit Shell Scripts**:
For Tiers 1-4, edit line 8 in tier\*.sh in each tier directory. Edit the output directory to your desired path. There must be a unique output directory for each tier.
 command="python ../../scripts/make-local-maps-IT.py.in -c ITpass2 -o /data/../OUTPUT\_DIRECTORY/tier\* $options"

For No Tier Cuts, edit line 8 in nocuts.sh in the nocuts directory. Edit the output directory to your desired path. There must be a unique output directory for each tier.
 command="python ../../scripts/make-local-maps-IT.py.in -c ITpass2 -o /data/../OUTPUT\_DIRECTORY/nocuts --ebins $options"

**Step 5. Run the Code**:
Navigate to 11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/run_all (or run_burnsample, for burnsample data). Make sure that the output directory is empty of any existing folders (delete all submit_20## files from previous runs), as the code does not overwrite any existing files, and instead skips them. Make sure you and in cobalt and have run cvmfs before running this script.

For tiers 1-4:

  cd tier#
  chmod +x tier#.sh
  ./tier#.sh

For nocuts:

  cd nocuts
  chmod +x nocuts.sh
  ./nocuts.sh

**Step 6. Iteration Method**:
*Note: this step may only be started once all files from the previous steps are done running and each directory has 360 files*

Combine fits files for each degree. The combine-local-maps.py.in script combines all the local fits files for each degree, producing 360 output files in a directory called 'combined'. These files are produced as described in the previous linked instructions, EXCEPT that the --reco flag is removed. Make sure you are in condor and have run cvmfs, then run the  following command is run in run_all/tier# or run_burnsample/tier#:

python ../../scripts/merge-reco.py.in -c ITpass2 -o /data/../OUTPUT\_DIRECTORY/tier# --submit_dir ./submit

The above line of code can take more than an hour to submit to condor, to run it in the background (so it continues to run if your computer turns off, e.g.), run 

nohup [command above] &

**Step 7. Reconstruction**:
*Note: this step may only be started once all the files from the previous step are done running and the combined directory has 360 files*

Enter Cobalt and run cvmfs 4.3.0

cd 11-year-Update-Paper-Figures/skymaps/cra-ahlers/build/configs

A .json config file must be created for each tier being reconstructed. You can copy the format needed by copying icecube.json from 11-year-Update-Paper-Figures/skymaps/cra-ahlers/build/examples into the config folder. For each copied .json file, make sure that thetamax is set to 55.0:

"thetamax": 55.0,

And "prefix" is changed to the output directory of your combined files:

"prefix": "/data/user/@USER_DIR@/@TIER@/ITpass2/combined/CR_ICECUBE_LOCAL_NSIDE64_degbin-",

Then run the following for each tier: ./bin/multi-llh --config ./configs/CONFIG_NAME.json  -o /data/@USER_DIR@/@TIER@/ITpass2/reco --iterations 20 --save-iter --smoothing-radius 20


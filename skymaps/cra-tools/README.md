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
On line 81, edit the submit directory to the path where your submit directory is:i
    default='/home/@USER@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/submit'
On line 223, update the cmd line to the location of your 'make-local-maps' file:
    cmd='/home/@USER@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/build/bin/make-local-maps'

The second is merge-reco.py.in. This script combines the local maps and uses the iteration method to estimate relative intensity and significance. Edit these lines:

On line 61, update the line to your output directory:
   default='/data/.../OUTPUT\_DIRECTORY',
On line 73, update the line to the version of cvmfs used in cmake:
   default='/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0',
On line 76, update the line to the path to your submit directory:
   default='/home/@USER@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/submit',
On line 214, update the line to the location of your 'combine-local-maps' folder:
   cmd  = '/home/@USER@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/build/bin/combine-local-maps'
On line 226, update the line to:
   cmd  = '/home/@USER@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/build/bin/illh-reco'
On line 228, update the line to:
   '--cfg','/home/@USER@/11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/build/resources/{detector}.json'.format(\*\*s\_opts),

We also have to edit 11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/private/SimpleDST.cc depending on which energy tier we intend to process:

For Tier 1 and No Tier Cuts (nocuts), comment out line 39 and uncomment line 40, so it reads
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
 command="python ../../scripts/make-local-maps-IT.py.in -c ITpass2 -o /data/../OUTPUT\_DIRECTORY/nocuts $options"

**Step 5. Run the Code**:
Navigate to 11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/run_all (or run_burnsample, for burnsample data). Make sure that the output directory is empty of any existing folders (delete all submit_20## files from previous runs), as the code does not overwrite any existing files, and instead skips them. Make sure you and in condor and have run cvmfs before running this script.

For tiers 1-4:

  cd tier#
  chmod +x tier#.sh
  ./tier#.sh

For nocuts:

  cd nocuts
  chmod +x nocuts.sh
  ./nocuts.sh

After you have processed daily fits files for each tier, you can use the script error_check.py to make sure no errors occurred in this process.
Navigate to 11-year-Update-Paper-Figures/skymaps/cra-tools/simpledst-maps/src/run_all/error_check.py (or .../src/run_burnsample/error_check.py).

On line 7 of error_checky.py, edit error_dir so it routes to 11-year-Update-Paper-Figures in your home directory.
 error_dir = Path('/home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps')

Once you have made this change run the script using the following command in terminal:
 python error_check.py

The contents of all error files (located in the submit directories in run_all/tier#/submit_####/npx4-error) will be printed. If there are no errors, nothing will be printed.
  
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

A .json config file must be created for each tier being reconstructed. You can copy the format needed by copying icetop.json from 11-year-Update-Paper-Figures/skymaps/cra-ahlers/build/examples into the config folder. For each copied .json file, make sure that thetamax is set to 55.0:

"thetamax": 55.0,

And "prefix" is changed to the output directory of your combined files:

"prefix": "/data/user/@USER_DIR@/@TIER@/ITpass2/combined/CR_ICECUBE_LOCAL_NSIDE64_degbin-",

Then run the following line in 11-year-Update-Paper/skymaps/cra-ahlers/build for each tier: ./bin/multi-llh --config ./configs/CONFIG_NAME.json  -o /data/@USER_DIR@/@TIER@/ITpass2/reco --iterations 20 --save-iter --smoothing-radius 20

**Step 8. Skymap Production**:
Navigate to the directory where your reconstruction fits files were produced. Change the name of the highest iteration of a file starting with CR_IceTop__64_360 to:
 combined_tX_iterationXX.fits.gz

And change the name of the highest iteration of a file starting with variance_IceTop__64_360 file to:
 significance_tX_iterationXX.fits.gz

Where X is replaced with the tier number and XX is replaced by the iteration number.

Do this for each tier. 
Copy each of these highest-iteration files into your 11-year-Update-Paper-Figures/skymaps/05_skymaps/output directory.

Navigate to 11-year-Update-Paper-Figures/skymaps/05_skymaps/2Dskymaps.ipynb . I recommend you access this script in JuypyterHub.

Edit the root directory to the location of your output files (in /data/ana/... or /data/user/...)
Edit the out directory to:
 home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps/05_skymaps/output
This is where the 2D skymap .png files will be output to.

Execute the blocks labeled #Relative Intensity Plots and #Significance Plots to create these maps for each energy tier.

**Step 9. 1D Projection Creation**
Navigate to 11-year-Update-Paper-Figures/skymaps/05_skymaps/alltiers.sh
Make sure each argument after $input corresponds to the name of the highest iteration of combined and significance reconstruction files. Make sure the input argument is linked to your renamed reconstruction files in 11-year-Update-Paper-Figures/skymaps/05_skymaps/output . Your output should be linked to the 1Ds folder in the same directory, this is where the 1 Dimensional Projection .png files will be produced. Save all changes in JupyterHub before you run this script in terminal so you are running the most updated version of the script. 
To run this script execute the following line in terminal:
 ./alltiers.sh
This will output .png files of the 1D projections of the 2D skymaps into your 1Ds folder.

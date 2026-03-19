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
cd ~/cra-tools-gb-edits/simpledst-maps
mkdir build
cd build
cmake ../src
make

**Step 2. Edit Code for the User**:
Navigate to cra-tools-gb-edits/simpledst-maps/src/scripts. We will edit 2 scripts here.

The first is make-local-maps-IT.py.in. This makes the local maps. Edit these lines:

On line 18, edit the fpath to the user directory where you have stored the burnsample root files.
   fpath = '/data/.../burnsamplemaps'
On line 63, edit the output directory to where you want the output directory to be:
   default='/data/.../OUTPUT\_DIRECTORY'
On line 75, edit the line to:
    default='/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1',
On line 78, edit the submit directory to the path where your submit directory is:
    default='/home/@USER_DIR@/cra-tools-gb-edits/simpledst-maps/src/submit'
On line 216, update the cmd line to the location of your 'make-local-maps' file:
    cmd='/home/@USER_DIR@/cra-tools-gb-edits/simpledst-maps/build/bin/make-local-maps'
The second is merge-reco.py.in. This script combines the local maps and uses the iteration method to estimate relative intensity and significance. Edit these lines:

On line 61, update the line to your output directory:
   default='/data/.../OUTPUT\_DIRECTORY',
On line 73, update the line to:
   default='/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1',
On line 76, update the line to the path to your submit directory:
   default='/home/USERNAME/cra-tools-gb-edits/simpledst-maps/src/submit',
On line 212, update the line to the location of your 'combine-local-maps' folder:
   cmd  = '/home/USERNAME/cra-tools-gb-edits/simpledst-maps/build/bin/combine-local-maps'
On line 224, update the line to:
   cmd  = '/home/USERNAME/cra-tools-gb-edits/simpledst-maps/build/bin/illh-reco'
On line 226, update the line to:
   '--cfg','/home/USERNAME/cra-tools-gb-edits/simpledst-maps/build/resources/{detector}.json'.format(\*\*s\_opts),
We also have to edit cra-tools-gb-edits/simpledst-maps/src/private/SimpleDST.cc depending on which energy tier we intend to process:

For Tier 1, comment out line 39 and uncomment line 40, so it reads
    reco = "ShowerPlane";
For Tiers 2, 3 or 4, uncomment out line 39, and comment line 40, so it reads
    reco = "Laputop";

**Step 3. Rebuild the Project**
If you have edited a .cc file, remake the project as follows:
cd ~/cra-tools-gb-edits/simpledst-maps
cd build
cmake ../src
make

**Step 4. Edit Shell Scripts**:
Edit line 8 in tier\*.sh in each tier directory. Edit the output directory to your desired path. There must be a unique output directory for each tier.
 command="python ../../scripts/make-local-maps-IT.py.in -c ITpass2 -o /data/../OUTPUT\_DIRECTORY/tier\* $options"

**Step 5. Run the Code**:
Navigate to cra-tools-gb-edits/simpledst-maps/src/run\_all. Make sure that the output directory is empty of any existing folders, as the code does not overwrite any existing files, and instead skips them.

  cd tier1
  chmod +x tier1.sh
  ./tier1.sh

**Step 6. Iteration Method**:
*Note: this step may only be started once all files from the previous steps are done running and each directory has 360 files*

Combine fits files for each degree. The combine-local-maps.py.in script combines all the local fits files for each degree, producing 360 output files in a directory called 'combined'. These files are produced as described in the previous linked instructions, EXCEPT that the --reco flag is removed, so the following command is run in run\_burnsample/tierX:
python ../../scripts/merge-reco.py.in -c ITpass2 -o /data/../OUTPUT\_DIRECTORY/tier1 --submit\_dir ./submit

**Step 7. Reconstruction**:
*Note: this step may only be started once all the files from the previous step are done running and the combined directory has 360 files*

Enter Cobalt and run cvmfs 4.3.0

cd cra-ahlers-llh-main/build/configs

Then run: ./bin/multi-llh --config ./configs/it\_t\*\_unblinded.json  -o /data/@USER\_DIR@/@TIER@/ITpass2/reco --iterations 20 --save-iter --smoothing-radius 20


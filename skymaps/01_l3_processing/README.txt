Instructions to process level3 files into .hdf5 and .root files for later analysis

For all steps, run in condor/submitter.
Step 1: Burnsample to HDF5
In 01_l3_processing, run CVMFS v4.3.0:

For the Burnsample (10% of the dataset), run:
python process_data_dag.py --root [OUTDIR] 

For the full dataset, run:
python process_data_dag.py -root [OUTDIR] --all_tiers

Change [OUTDIR] to the address in \data\user\@USER_DIR@ where you want your root files to be output to 

Wait until the job is done by checking the condor queue (condor_q) before moving on to the next step

To check if any errors were produced in this process:

Navigate to 11-year-Update-Paper-Figures/skymaps/01_l3_processing/error_check_root.py. Edit the error_dir path in line 7 to  your user directory.
 error_dir = Path('/home/@USER_DIR@/11-year-Update-Paper-Figures/skymaps')

Run this script in terminal using the command:
 python error_check_root.py

This script will print the the error files and their contents that have thrown errors. If no errors occurred, nothing will be printed. 
 
 

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


 

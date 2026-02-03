Instructions to process level3 files into .hdf5 and .root files for later analysis

For all steps, run in condor/submitter.
Step 1: Burnsample to HDF5
In 01_l3_processing, run CVMFS v4.3.0:

Then run:
python process_data_dag.py --root [OUTDIR] --all_data ALL_DATA 

Wait until the job is done by checking the condor queue (condor_q) before moving on to the next step


 

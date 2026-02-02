Instructions to process level3 files into .hdf5 and .root files for later analysis

For all steps, run in condor/submitter.
Step 1: Burnsample to HDF5
In 01_l3_processing, run CVMFS v4.3.0:

Then run:
python process_data_dag.py --hdf [output dir] 

Step 2: Burnsample to ROOT
In 01_l3_processing, run CVMFS v4.3.0:

Then run:
python process_data_dag.py --root [output dir] 

The two steps can be combined (e.g. process_data_dag.py --hdf [dir1] --root [dir2]) or run for specific years (e.g. process_data_dag.py -y 2011 2013).



 

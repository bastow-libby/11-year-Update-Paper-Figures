Instructions to process simulation files into .hdf5 files for later analysis

For all steps, run in condor/submitter.

Step 1: Simulation to HDF5
In sim_hdf5, run CVMFS v4.3.0:

Then run:
python run_stuff.py -o [output dir] -l [log dir] -y [years]

Step 2: Enter the sim_study directory and follow the instructions in the README file to create the simulation arrays and plots 

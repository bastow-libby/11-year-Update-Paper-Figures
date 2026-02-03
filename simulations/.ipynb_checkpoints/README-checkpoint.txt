Instructions to process simulation files into .hdf5 files for later analysis

For all steps, run in condor/submitter.
Make sure the line "combo() { ${CVMFS}/py3-v4.3.0/icetray-env /icetray/v1.9.2; }" is in your source file

Step 1: Simulation to HDF5

In sim_hdf5, run CVMFS v4.3.0:

Then run:
python run_stuff.py -o [output dir] -l [log dir] -y [years]

Step 2: Enter the sim_study directory and follow the instructions in the README file to create the simulation arrays and plots 

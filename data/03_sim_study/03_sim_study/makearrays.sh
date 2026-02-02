#!/bin/bash
dir=/data/user/rchapagain/cra_analysis/
#python sw1_nstation_binning.py $dir"sim_2018" -o ./arrays/2018/default -e [[6,13]]
#python sw1_nstation_binning.py $dir"sim_2018"  -o ./arrays/2018/1by3 -e [[6,13]]
#python sw1_nstation_binning.py $dir"sim_2018" -o ./arrays/2018/2by3 -e [[7,14]]


#python sw1_nstation_binning.py $dir"sim_2015" -o ./arrays/2015/ns
#python sw1_nstation_binning.py $dir"sim_2018" -o ./arrays/2018/ns

python sw1_nstation_binning.py $dir"sim_2012_simweights" -o ./arrays/simweight

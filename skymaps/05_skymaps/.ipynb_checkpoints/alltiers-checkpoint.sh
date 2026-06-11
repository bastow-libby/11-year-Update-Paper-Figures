#!/bin/bash

input='./anti_output/'
output='./1Ds_anti/'
fit='--gaussian'
fit2='--fit 3'
fit3='--fit 1'
fit4='--flat'
range3='-M 3 -m -3'
range2='-M 1 -m -1'
range1='-M 1.5 -m -1.5'
range4='-M 4 -m -4'
range5='-M 2.5 -m -2.5'
range6='-M 3.5 -m -3.5'
range7='-M 2 -m -2'
python 1d_proj.py $input"combined_t1_iteration20.fits.gz" --relerr $input"significance_t1_iteration20.fits.gz" -o $output"1d_t1_anti.png" $fit3 $fit4 $range2 -d -90 -D -40 -z -S 3 --title "Tier 1: 300 TeV 2011-14"
#"No Cuts Sidereal: 2011-21"#
#"Tier 1 Solar: 310 TeV 2011-14"
python 1d_proj.py $input"combined_t2_iteration03.fits.gz" --relerr $input"significance_t2_iteration03.fits.gz" -o $output"1d_t2_anti.png" $fit3 $fit4 $range1 -d -90 -D -40 -z -S 3 --title "Tier 2: 900 TeV 2011-14"
python 1d_proj.py $input"combined_t3_iteration20.fits.gz" --relerr $input"significance_t3_iteration20.fits.gz" -o  $output"1d_t3_anti.png" $fit3 $fit4 $range1 -d -90 -D -40 -z -S 3 --title "Tier 3: 2.3 PeV 2011-21"
python 1d_proj.py $input"combined_t4_iteration20.fits.gz" --relerr $input"significance_t4_iteration20.fits.gz" -o $output"1d_t4_anti.png" $fit3 $fit4 $range4 -d -90 -D -40 -z -S 3 --title "Tier 4: 6.6 PeV 2011-21"
#Tier 4: 2011-21 6.6 PeV
#python 1d_proj.py $input"combined_nc_iteration20.fits.gz" --relerr $input"significance_nc_iteration20.fits.gz" -o $output"1d_nc.png" $fit2 $range1 -d -90 -D -36 -z -S 3 --title "No Tier Cuts 2011-21"
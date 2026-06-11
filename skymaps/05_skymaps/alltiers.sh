#!/bin/bash

input='./output'
output='./1Ds/sidereal_unblinded/'

fit='--gaussian'
fit1='--fit 1'
fit3='--fit 3'
fit4='--flat'

range1='-M 1 -m -1'
range2='-M 1.5 -m -1.5'
range3='-M 2 -m -2'
range4='-M 2.5 -m -2.5'
range5='-M 3 -m -3'
range6='-M 3.5 -m -3.5'
range7='-M 4 -m -4'

#tier 1 
python 1d_proj.py $input"combined_t1_iteration20.fits.gz" --relerr $input"significance_t1_iteration20.fits.gz" -o $output"1d_t1.png" $fit2 $range2 -d -90 -D -40 -z -S 3 --title "Tier 1: 300 TeV 2011-14"

#tier 2 
python 1d_proj.py $input"combined_t2_iteration20.fits.gz" --relerr $input"significance_t2_iteration20.fits.gz" -o $output"1d_t2.png" $fit2 $range3 -d -90 -D -40 -z -S 3 --title "Tier 2: 900 TeV 2011-14"

#tier 3 
python 1d_proj.py $input"combined_t3_iteration20.fits.gz" --relerr $input"significance_t3_iteration20.fits.gz" -o  $output"1d_t3.png" $fit3 $fit4 $range1 -d -90 -D -40 -z -S 3 --title "Tier 3: 2.3 PeV 2011-21"

#tier 4 
python 1d_proj.py $input"combined_t4_iteration20.fits.gz" --relerr $input"significance_t4_iteration20.fits.gz" -o $output"1d_t4.png" $fit3 $fit4 $range4 -d -90 -D -40 -z -S 3 --title "Tier 4: 6.6 PeV 2011-21"

#nocuts 
python 1d_proj.py $input"combined_nc_iteration20.fits.gz" --relerr $input"significance_nc_iteration20.fits.gz" -o $output"1d_nc.png" $fit2 $range1 -d -90 -D -36 -z -S 3 --title "No Tier Cuts 2011-21"

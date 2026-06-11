#!/usr/bin/env python

import json
from pathlib import Path

if __name__ == "__main__":
    error_dir = Path('/home/@USER@/11-year-Update-Paper-Figures/skymaps')

    temp = error_dir / f'01_l3_processing/submit_root'
    files = sorted(temp.glob('20??/npx4-error/*.error'))
    #rootprod_2017051_25376.error'))
    
    truelength = 98
    #print(truelength)

    for file_path in files:

        with open(file_path, 'r') as f:
            lines = f.readlines()
            corr_line = ['NOTICE (I3Tray): I3Tray finishing... (I3Tray.cxx:526 in void I3Tray::Execute(bool, unsigned int))\n']

            if str(lines) != str(corr_line): 
                print(str(file_path))
                print("                                              ")
                for line in lines:
                    print(f' {line}')
                print("----------------------------------------------")
                print("                                              ")


                

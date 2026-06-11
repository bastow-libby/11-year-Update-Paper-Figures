#!/usr/bin/env python

from pathlib import Path

if __name__ == "__main__":

    error_dir = Path('/home/amcclure/11-year-Update-Paper-Figures/skymaps')

    tiers = ["tier1", "tier2", "tier3", "tier4", "nocuts"]
    
    for tier in tiers:

        temp = error_dir / f'cra-tools/simpledst-maps/src/run_burnsample/{tier}'
        files = sorted(temp.glob('submit_????/npx4-error/*.error'))

        for file_path in files:

            with open(file_path, 'r') as f:
                lines = f.readlines()

            if lines != []:
                print(str(file_path))
                for line in lines:
                    print(f'  {line}')

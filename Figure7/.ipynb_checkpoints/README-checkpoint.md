## Folder for creating Figure 7 (Angular Power Spectrum) in the 11 years with IceTop paper
### Scripts in this Folder
- plotMaker.py : wrapper script that can generate both the uncertainties and the angular power spectrum with simplified commands.
  - How to run: python [code] -f [input file path] -t t# -o [output file path] (recommended first run: -m to make the uncertainties.)
  - input and output file paths have set defualts, change these for your needs.
- aps.py : generates and plots the angular power spectrum
  - How to run: python [code] -f [input file path] -i [isoErr path] --staterr [statErr path] --syserr [sysErr path] -o [output file path] -l [label]
- isoErr.py : generates the isotropic bands
  - How to run: python [code] -f [input file path] -o [output file path] (optional: -n [int amount of times to run])
- statErr.py : generates the statistical error bars
  - How to run: python [code] -f [input file path] -o [output file path] (optiona: -n [int amount of times to run])
- sysErr.py : generates the systematic error bars
  - How to run: python [code] -f [input file path] -o [output file path] (optional: -n [int amount of times to run])
- map_functions.py: contains functions used in the above files for plot creation.

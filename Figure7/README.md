## Folder for creating Figure 7 (Angular Power Spectrum) in the 11-year IceTop anisotropy paper

### How to Run
- if you do not have access to existing uncertainty files, produce them using
  submit_err.py
  - submit_err.py will output a series of fake Cl values. Use calc_err.py to
    simplify these to the uncertainties that the plotting script looks for
- run plot_maker.py to replicate paper figures

### Scripts in this Folder
- aps.py : generates and plots the angular power spectrum
- calc_err.py : combines fake Cl values from stat, sys, and isoErr.py to
  calculate uncertainty values
- isoErr.py : generates fake Cl values for isotropic noise bands
- plot_maker.py : wrapper script that generates the power spectra images
- README.md : this file
- statErr.py : generates fake Cl values for statistical error bars
- submit_err.py : wrapper script for cluster submission of iso/stat/sysErr
- submitter/ : directory with tools for cluster submission
- sysErr.py : generates fake Cl values for systematic error bars

#Shell script to generate all angular power spectra in their respective folders
python3 plotMaker.py -f /data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output/outpute/finalcombinedfits -t 1 -o . -m -s 0 -l 300_TeV -il
python3 plotMaker.py -f /data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output/outpute/finalcombinedfits -t 2 -o . -m -s 0 -l 900_TeV
python3 plotMaker.py -f /data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output/outpute/finalcombinedfits -t 3 -o . -m -s 0 -l 2300_TeV
python3 plotMaker.py -f /data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output/outpute/finalcombinedfits -t 4 -o . -m -s 0 -l 6600_TeV

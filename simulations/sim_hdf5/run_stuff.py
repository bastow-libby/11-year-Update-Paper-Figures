import os, fileinput, sys
import argparse

fnum_dict = {2012: ("SIBYLL2.1", [12360, 12630, 12631, 12362]),
             2015: ("SIBYLL2.1", [20174, 20178, 20179, 20180]),
             2018: ("SIBYLL2.3d", [22570, 22580, 22583, 22586])}

pars = ["p", "He", "O", "Fe"]

def runSim(outfile, filepath):
    for line in fileinput.input('steeringcard_dstsim',inplace=1):
	    arg_string = outfile + ' '  + filepath
	    if 'arguments' in line:
	        line = line.replace(line, 'arguments = '+arg_string+'\n')
	    sys.stdout.write(line)
    os.system('condor_submit steeringcard_dstsim')
    print('running for', filepath)


# Main method
if __name__ == "__main__":
    # Set up command line options
    parser = argparse.ArgumentParser(description=
                                     "Process sim l3.i3 files as hdf5 files.")
    parser.add_argument("-o", "--output", dest="output", type=str, default=None,
                        help="Output directory for hdf5 files")
    parser.add_argument("-l", "--log", dest="log", type=str,
                        default=None, help="Directory for condor output/error/log files")
    parser.add_argument("-y", "--years", dest="years", type=int, nargs="+", choices=[2012, 2015, 2018],
                        default=None, help="Simulation year(s) to process out of 2012, 2015, and/or 2018")
    parser.add_argument("-p", "--particles", dest="particles", nargs="+", type=str,
                        choices= ["p", "He", "O", "Fe"],
                        default= ["p", "He", "O", "Fe"],
                        help="Particles to process (default: p He O Fe)")
    args = parser.parse_args()
    
    # Create directories where necessary
    if not os.path.isdir(args.output):
        print('Creating output directory {}'.format(args.output))
        os.makedirs(args.output)
    if not os.path.isdir(args.log):
        print('Creating condor output directory {}'.format(args.log))
        os.makedirs(args.log)

    for line in fileinput.input('steeringcard_dstsim',inplace=1):
        if 'output =' in line:
            line = line.replace(line,'output = '+ args.log +'/condor.out.$(ClusterId)\n')
        elif 'error =' in line:
            line = line.replace(line,'error = '+ args.log +'/condor.err.$(ClusterId)\n')
        elif 'log =' in line:
            line = line.replace(line,'log = '+ args.log +'/condor.log.$(ClusterId)\n')
        sys.stdout.write(line)

    for year in args.years:
        model, fnums = fnum_dict[year]
        path = "/data/ana/CosmicRay/IceTop_level3/sim/IC86.{year}/{model}/".format(year = year, model = model)
        file = "Level3_IC86.{year}_{model}_*.i3.*".format(year = year, model = model)
        for particle in args.particles:
            ind = pars.index(particle)
            if year == 2018:
                mid_path = "{particle}_allE_links/".format(particle = particle)
            else:
                mid_path = "{particle}/{fnum}_v1s/".format(particle = particle, fnum = fnums[ind])
            path_full = path + mid_path + file 
            outpath = '{outdir}/l3_{year}_{part}.hdf'.format(outdir = args.output, year = year, part = particle)
            runSim(outpath, path_full)

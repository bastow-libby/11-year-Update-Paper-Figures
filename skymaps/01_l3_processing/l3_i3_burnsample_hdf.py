#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT: /icetray/v1.9.2
from icecube import icetray, dataclasses, simclasses, dataio, tableio, toprec

from icecube.hdfwriter import I3HDFWriter

from icecube.icetray import I3Tray
from icecube.icetop_Level3_scripts.functions import count_stations

import sys, os
import numpy as np
from argparse import ArgumentParser
from glob import glob

def add_IceTop_quality_cuts(frame):
        passed = all(frame['IT73AnalysisIceTopQualityCuts'].values())
        #print "test " + str(passed)
        frame['passed_IceTopQualityCuts'] = icetray.I3Bool(passed)


def add_nstations(frame, pulses='IceTopHLCSeedRTPulses_SnowUnAttenuated'):
    nstation = count_stations(dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulses))
    frame['NStations'] = icetray.I3Int(nstation)


parser = ArgumentParser(description=__doc__)
parser.add_argument('outdir', type=str)
parser.add_argument('run_year', type=str)
parser.add_argument('year', type=str)
parser.add_argument('month', type=str)
#parser.add_argument('infile', nargs='*')
#parser.add_argument('outfile')
opts = parser.parse_args()

'''
run_year = int(sys.argv[1])
year = int(sys.argv[2])
month = int(sys.argv[3])

#physics = glob("/data/exp/IceCube/"+year+"/filtered/level2/*/Level2_IC86_corsika_icetop.010410.*.i3.bz2")
'''
physics = glob("/data/ana/CosmicRay/IceTop_level3/exp/IC86.{run_year}_pass2_v0*/{cal_year}/{month}*/*/Level3_IC86.{run_year}*0_Subrun*.i3.*".format(run_year = opts.run_year,
                                                                                                                                                   cal_year = opts.year,
                                                                                                                                                   month = opts.month))
#notice that this only takes events from runs ending in 0 to only gather 10% of events

#physics = glob("/data/ana/CosmicRay/IceTop_level3/exp/IC86."+str(run_year)+"/"+str(year)+"/*/Level3_IC86."+str(run_year)+"*0_Subrun*.i3.*")

#/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/"+dataset+"/Level3_IC86.2012_"+dataset+"_Run0*.i3.gz")    # glob() the list of files from the disk

physics.sort()
 
#for i in physics:
#	print(i)
if len(physics)==0:
	print("This didn't work - no files globbed")

#gcd = '/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/GCD/Level3_'+dataset+'_GCD.i3.gz'
#gcd = '/home/mplum/GeoCalibDetectorStatus_2012.56063_V1_OctSnow_scint.i3'
subevent = 'IceTopSplit'

temp_outfile = opts.outdir+'l3_data_run_config_{run_year}_year_{cal_year}_{month}.hdf5'.format(run_year = opts.run_year, cal_year = opts.year, month = opts.month)
if len(physics)!=0:
	tray = I3Tray()
	#tray.AddModule('I3Reader', 'reader', filenamelist=opts.infile)
	tray.AddModule('I3Reader', 'reader', DropBuffers= False, filenamelist=physics)

	tray.Add(add_nstations,
		 If=lambda frame: 'IceTopHLCSeedRTPulses_SnowUnAttenuated' in frame)

	tray.Add(add_IceTop_quality_cuts,
		 If=lambda frame: 'IT73AnalysisIceTopQualityCuts' in frame)


	tray.AddSegment(I3HDFWriter, 'scribe',
	    #Output=opts.outfile,
	    Output=temp_outfile,
	    Keys=[
		  'I3EventHeader'
		  'QFilterMask',
		  'Laputop',
		  'LaputopParams',
		  'ShowerPlane',
		  'ShowerPlaneParams',
		  'IT73AnalysisIceTopQualityCuts',
		  'NStations'
		  ],
	    Types=[],
	    SubEventStreams=[subevent],
	    #DropOrphanStreams=[icetray.I3Frame.DAQ],
	    #streams=[icetray.I3Frame.Physics],
	)


	tray.Execute()																															


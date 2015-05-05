#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 3:
	print "This script must be called with three arguments!"
	print "Usage:\n   palmaprocesslauncher.py file time"
	sys.exit(1)

pattern = re.search('.*_([0-9]+)\.cfg',sys.argv[1])
if pattern != None:
	scriptfile = open("fermions"+pattern.group(1)+".cmd",'w')
	scriptdata = "#PBS -o output.dat \n\
#PBS -l walltime="+sys.argv[2]+":00:00,nodes=1:westmere:ppn=12 \n\
#PBS -A p0muenst \n\
#PBS -M spiem_01@uni-muenster.de \n\
#PBS -m ae \n\
#PBS -q default \n\
#PBS -N fermion"+pattern.group(1)+" \n\
#PBS -j oe \n\
\n\
cd $PBS_O_WORKDIR \n\
\n\
./leonardQCD.exe --configfile="+sys.argv[1]+" > logfile"+pattern.group(1)+".txt\n"
	scriptfile.write(scriptdata);
	

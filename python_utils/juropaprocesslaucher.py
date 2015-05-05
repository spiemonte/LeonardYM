#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 3:
	print "This script must be called with three arguments!"
	print "Usage:\n   juropaprocesslauncher.py file time"
	sys.exit(1)

pattern = re.search('.*_([0-9]+)\.cfg',sys.argv[1])
if pattern != None:
	scriptfile = open("fermions"+pattern.group(1)+".cmd",'w')
	scriptdata = "#!/bin/bash -x \n\
#MSUB -l nodes=8:ppn=8 \n\
#MSUB -e /lustre/jhome3/hhh04/hhh041/WilsonFlow/err_first.txt \n\
#MSUB -o /lustre/jhome3/hhh04/hhh041/WilsonFlow/out_first.txt \n\
#MSUB -v tpt=8 \n\
#MSUB -l walltime="+sys.argv[2]+":00:00 \n\
### start of jobscript \n\
export OMP_NUM_THREADS=8 \n\
cd $PBS_O_WORKDIR \n\
echo \"workdir: $PBS_O_WORKDIR\" \n\
 \n\
# NSLOTS = nodes * ppn / tpt = 8 * 8 / 8 = 8 \n\
NSLOTS=8 \n\
mpiexec -np $NSLOTS --exports=OMP_NUM_THREADS ./leonardQCD.exe --configfile="+sys.argv[1]+" > logfile"+pattern.group(1)+".txt \n"
	scriptfile.write(scriptdata);
	

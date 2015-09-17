#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 4 and len(sys.argv) != 5:
	print "This script must be called with four/five arguments!"
	print "Usage:\n   palmaprocesslauncher.py file time nodes"
	print "Usage:\n   palmaprocesslauncher.py file time nodes load_layout"
	sys.exit(1)

load_layout = ""

if len(sys.argv) == 5:
	load_layout = "--load_layout"
	
pattern = re.search('.*_([0-9]+)\.cfg',sys.argv[1])
if pattern != None:
	scriptfile = open("script"+pattern.group(1)+".cmd",'w')
	scriptdata = "#!/bin/bash \n\
# DO NOT USE environment = COPY_ALL \n\
#@ job_type = parallel \n\
#@ class = fat \n\
#@ node = "+sys.argv[3]+" \n\
#@ tasks_per_node = 4 \n\
#@ wall_clock_limit = "+sys.argv[2]+":00:00 \n\
#@ job_name = fermion"+pattern.group(1)+" \n\
#@ network.MPI = sn_all,not_shared,us \n\
#@ initialdir = "+os.getcwd()+" \n\
#@ output = job$(jobid).out \n\
#@ error = job$(jobid).err \n\
#@ notification=error \n\
#@ notify_user=spiemonte@uni-muenster.de \n\
#@ queue \n\
. /etc/profile \n\
. /etc/profile.d/modules.sh \n\
#setup of environment \n\
export MP_SINGLE_THREAD=no \n\
export OMP_NUM_THREADS=10 \n\
# Pinning \n\
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS \n\
#optional:  \n\
#module load mpi_pinning/hybrid_blocked \n\
module load gcc/4.9 \n\
mpiexec -n "+str(int(sys.argv[3])*4)+" ./leonardQCD.exe --configfile="+sys.argv[1]+" "+load_layout+" > logfile"+pattern.group(1)+".txt \n\
 \n\
echo -e \"configfile = "+sys.argv[1].replace("./","")+" \\nscriptfile = script"+pattern.group(1)+".cmd\\n\" > script"+pattern.group(1)+"_reloadme.txt \n"
	scriptfile.write(scriptdata);



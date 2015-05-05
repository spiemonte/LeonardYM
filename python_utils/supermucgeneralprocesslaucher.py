#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 4:
	print "This script must be called with four arguments!"
	print "Usage:\n   palmaprocesslauncher.py file time nodes"
	sys.exit(1)

pattern = re.search('.*_([0-9]+)\.cfg',sys.argv[1])
if pattern != None:
	scriptfile = open("fermions"+pattern.group(1)+".cmd",'w')
	scriptdata = "#!/bin/bash \n\
## \n\
#@ energy_policy_tag = spiemonte \n\
#@ max_perf_decrease_allowed = -10 \n\
## \n\
#@ wall_clock_limit = "+sys.argv[2]+":00:00 \n\
#@ job_name = fermion"+pattern.group(1)+" \n\
#@ job_type = parallel \n\
#@ class = general \n\
#@ node = "+sys.argv[3]+" \n\
#@ total_tasks = "+sys.argv[3]+" \n\
#@ island_count = 1 \n\
#@ node_usage = not_shared \n\
#@ network.MPI = sn_all,not_shared,us \n\
#@ initialdir = "+os.getcwd()+" \n\
#@ output = job$(jobid).out \n\
#@ error = job$(jobid).err \n\
#@ notification=always \n\
#@ notify_user=spiemonte@uni-muenster.de \n\
#@ queue \n\
\n\
. /etc/profile \n\
. /etc/profile.d/modules.sh \n\
\n\
module unload mpi.ibm \n\
module load mpi.ibm/1.3_gcc \n\
\n\
export MP_SINGLE_THREAD=no \n\
export OMP_NUM_THREADS=32 \n\
# Pinning \n\
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS \n\
\n\
mpiexec -n "+sys.argv[3]+" ./leonardQCD.exe --configfile="+sys.argv[1]+" > logfile"+pattern.group(1)+".txt\n"
	scriptfile.write(scriptdata);
	

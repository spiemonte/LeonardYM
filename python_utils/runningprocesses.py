#!/usr/bin/python

import commands
import re
import sys
import os


if len(sys.argv) != 2:
	print "This script must be called with two arguments!"
	print "Usage:\n   runningprocesses.py user"
	sys.exit(1)
	
rawdirs = commands.getstatusoutput("qstat -f1 | grep "+sys.argv[1]+" | grep init_work_dir")[1].split('\n')
rawfiles = commands.getstatusoutput("qstat | grep "+sys.argv[1])[1].split('\n')

dirs = []
files = []

for line in rawdirs:
	pattern = re.search("init_work_dir = ([\/A-Za-z0-9\_]+)",line)
	dirs.append(pattern.group(1))
	
for line in rawfiles:
	pattern = re.search("([0-9A-Za-z]+\.cmd)", line)
	files.append(pattern.group(1))

for namefile in os.listdir(os.getcwd()):
	pattern = re.search('([0-9A-Za-z]+\.cmd)',namefile)
	if pattern != None:
		flag = True
		for i in range(0,len(files)-1):
			if os.getcwd()+"/"+pattern.group(1) == str(dirs[i]) + "/" + str(files[i]):
				flag = False
		if flag:
			print "File not running " + os.getcwd()+"/"+pattern.group(1)

for i in range(0,len(files)-1):
	print "Running file " + str(dirs[i]) + "/" + str(files[i]) + " " + rawfiles[i].split()[4]
	


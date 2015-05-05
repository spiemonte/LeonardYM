#!/usr/bin/python

import re
import os
import sys
import time
from subprocess import call

if len(sys.argv) != 2:
	print "This script must be called with three arguments!"
	print "Usage:\n   reloadprocesses.py folder"
	sys.exit(1)
	
while True:
	for namefile in os.listdir(sys.argv[1]):
		pattern = re.search('[a-zA-Z0-9\-\_]\_reloadme.txt',namefile)
		if pattern != None:
			inputfile = open(sys.argv[1]+namefile,'r')
			filedata = inputfile.read()
			print filedata
			patternconfig = re.search('configfile\ *=\ *([a-zA-Z0-9\_\.]+)',filedata)
			patternscript = re.search('scriptfile\ *=\ *([a-zA-Z0-9\_\.]+)',filedata)
			call(["./restartsimulation.py ", sys.argv[1]+patternconfig.group(1)])
			time.sleep(5)
			call(["msub ", sys.argv[1]+patternscript.group(1)])
			print  sys.argv[1]+namefile
			os.remove(sys.argv[1]+namefile)
	time.sleep(10)

#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 3:
	print "This script must be called with two arguments!"
	print "Usage:\n   changekappa.py file new_number_warmup_sweeps"
	sys.exit(1)
	
inputfile = open(sys.argv[1],'r')

filedata = inputfile.read()
pattern = re.search('(number_warm_up_sweeps.*=.*[0-9]+)',filedata)

if pattern != None:
	filedata = filedata.replace(pattern.group(1),'number_warm_up_sweeps = '+sys.argv[2])
else:
	print "Warning, no number_warm_up_sweeps found!"

inputfile.close()

outputfile = open(sys.argv[1],'w')
outputfile.write(filedata)
outputfile.close()





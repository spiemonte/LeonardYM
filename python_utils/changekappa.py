#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 3:
	print "This script must be called with two arguments!"
	print "Usage:\n   changekappa.py file newkappa"
	sys.exit(1)
	
newkappa = float(sys.argv[2])
inputfile = open(sys.argv[1],'r')

filedata = inputfile.read()
pattern = re.search('kappa\ *=\ *([0-9]+\.[0-9]+)',filedata)

if pattern != None:
	replacerule = re.search('(kappa\ *=\ *[0-9]+\.[0-9]+)',filedata)
	filedata = filedata.replace(replacerule.group(1),'kappa = '+sys.argv[2])
	filedata = filedata.replace(str(float(pattern.group(1))).replace("0.","")+"k",str(float(newkappa)).replace("0.","")+"k")
else:
	print "Warning, no kappa found!"

inputfile.close()

outputfile = open(sys.argv[1],'w')
outputfile.write(filedata)
outputfile.close()





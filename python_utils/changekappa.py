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
pattern = re.search('kappa.*=.*([0-9]+\.[0-9]+)',filedata)

if pattern != None:
	replacerule = re.search('(kappa.*=.*[0-9]+\.[0-9]+)',filedata)
	filedata = filedata.replace(replacerule.group(1),'kappa = '+sys.argv[2])
	oldf = 'k'
	newf = 'k'
	if float(pattern.group(1)) < 0.1:
		oldf = '0' + str(int(float(pattern.group(1))*10000)) + oldf
	else:
		oldf = str(int(float(pattern.group(1))*10000)) + oldf
	if newkappa < 0.1:
		newf = '0' + str(int(newkappa*10000)) + newf
	else:
		newf = str(int(newkappa*10000)) + newf
	filedata = filedata.replace(oldf,newf)
else:
	print "Warning, no kappa found!"

inputfile.close()

outputfile = open(sys.argv[1],'w')
outputfile.write(filedata)
outputfile.close()





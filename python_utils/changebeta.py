#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 3:
	print "This script must be called with two arguments!"
	print "Usage:\n   changebeta.py file newbeta"
	sys.exit(1)
	
newbeta = float(sys.argv[2])
inputfile = open(sys.argv[1],'r')

filedata = inputfile.read()
pattern = re.search('beta\ *=\ *([0-9]+\.[0-9]+)',filedata)

if pattern != None:
	replacerule = re.search('(beta\ *=\ *[0-9]+\.[0-9]+)',filedata)
	filedata = filedata.replace(replacerule.group(1),'beta = '+sys.argv[2])
	oldf = str(int(float(pattern.group(1))*1000))+'b'
	newf = str(int(newbeta*1000))+'b'
	filedata = filedata.replace(oldf,newf)
else:
	print "Warning, no beta found!"

inputfile.close()

outputfile = open(sys.argv[1],'w')
outputfile.write(filedata)
outputfile.close()





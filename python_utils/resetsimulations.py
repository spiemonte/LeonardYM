#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 2:
	print "This script must be called with two arguments!"
	print "Usage:\n   resetsimulations.py file"
	sys.exit(1)
	
inputfile = open(sys.argv[1],'r')

filedata = inputfile.read()

patternoutputfile = re.search('output_name\ *=\ *(.+_)([0-9]+)',filedata)

patternjump = re.search('Output,\ *([0-9]+)',filedata)

patternstart = re.search('(start.*=.*readstart)',filedata)

if patternstart != None:
	filedata = filedata.replace(patternstart.group(1),'start = hotstart')
	
numberLastConf = 0

patterninputnumber = re.search('(input_number.*=\ *.+)',filedata)
if patterninputnumber != None:
	filedata = filedata.replace(patterninputnumber.group(1),'input_number = 0')

patternoutputnumber = re.search('(output_offset.*=\ *.+)',filedata)
if patternoutputnumber != None:
	if patternjump != None:
		numberLastConf += int(patternjump.group(1))
	filedata = filedata.replace(patternoutputnumber.group(1),'output_offset = 0')
	
if patternoutputfile != None:
	newnumber = 1
	filedata = filedata.replace(patternoutputfile.group(1)+patternoutputfile.group(2),patternoutputfile.group(1)+str(newnumber))

inputfile.close()

outputfile = open(sys.argv[1],'w')
outputfile.write(filedata)
outputfile.close()





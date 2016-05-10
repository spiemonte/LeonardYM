#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 4:
	print "This script must be called with three arguments!"
	print "Usage:\n   changeiodirectory.py file newdirectoryconfigs newdirectorymeas"
	sys.exit(1)
	
inputfile = open(sys.argv[1],'r')

filedata = inputfile.read()

patterninput = re.search('input_directory_configurations.*=.*(/.+)',filedata)
patternoutput = re.search('output_directory_configurations.*=.*(/.+)',filedata)

patternmeas = re.search('output_directory_measurements.*=.*(/.+)',filedata)

if patterninput != None:
	replacerule = re.search('(input_directory_configurations.*=.*/.+)',filedata)
	filedata = filedata.replace(replacerule.group(1),'input_directory_configurations = '+sys.argv[2])

if patternoutput != None:
	replacerule = re.search('(output_directory_configurations.*=.*/.+)',filedata)
	filedata = filedata.replace(replacerule.group(1),'output_directory_configurations = '+sys.argv[2])
	
if patternmeas != None:
	replacerule = re.search('(output_directory_measurements.*=.*/.+)',filedata)
	filedata = filedata.replace(replacerule.group(1),'output_directory_measurements = '+sys.argv[3])

inputfile.close()

outputfile = open(sys.argv[1],'w')
outputfile.write(filedata)
outputfile.close()





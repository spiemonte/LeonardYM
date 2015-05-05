#!/usr/bin/python

import re
import os
import sys

if len(sys.argv) != 2:
	print "This script must be called with two arguments!"
	print "Usage:\n   restartsimulations.py file"
	sys.exit(1)
	
inputfile = open(sys.argv[1],'r')

filedata = inputfile.read()

patterninput = re.search('input_directory_configurations\ *=\ *([a-zA-Z\_\/0-9\-]+)',filedata)
patternoutput = re.search('output_directory_configurations\ *=\ *([a-zA-Z\_\/0-9\-]+)',filedata)
patternoutputname = re.search('output_configuration_name\ *=\ *([a-zA-Z\_\/0-9\-]+)',filedata)
patternoutputfile = re.search('output_name\ *=\ *(.+_)([0-9]+)',filedata)

patternjump = re.search('Output,\ *([0-9]+)',filedata)

patternstart = re.search('(start.*=.*hotstart)',filedata)

if patternstart != None:
	filedata = filedata.replace(patternstart.group(1),'start = readstart')
	
numberLastConf = 0

print patternoutput.group(1)

for namefile in os.listdir(patternoutput.group(1)):
	pattern = re.search('0:([0-9]+):'+patternoutputname.group(1),namefile)
	if pattern != None:
		tmp = int(pattern.group(1))
		if tmp > numberLastConf:
			numberLastConf = tmp

if numberLastConf == 0:
	print "Warning, no old configurations found, trying with leonard format"
	for namefile in os.listdir(patternoutput.group(1)):
		pattern = re.search(patternoutputname.group(1)+'\_([0-9]+)',namefile)
		if pattern != None:
			tmp = int(pattern.group(1))
			if tmp > numberLastConf:
				numberLastConf = tmp

print "Restarting from configuration number: " + str(numberLastConf)

patterninputnumber = re.search('(input_number.*=\ *.+)',filedata)
if patterninputnumber != None:
	filedata = filedata.replace(patterninputnumber.group(1),'input_number = '+str(numberLastConf))

patternoutputnumber = re.search('(output_offset.*=\ *.+)',filedata)
if patternoutputnumber != None:
	if patternjump != None:
		numberLastConf += int(patternjump.group(1))
	filedata = filedata.replace(patternoutputnumber.group(1),'output_offset = '+str(numberLastConf))
	
if patternoutputfile != None:
	newnumber = int(patternoutputfile.group(2))+1
	filedata = filedata.replace(patternoutputfile.group(1)+patternoutputfile.group(2),patternoutputfile.group(1)+str(newnumber))

inputfile.close()

outputfile = open(sys.argv[1],'w')
outputfile.write(filedata)
outputfile.close()





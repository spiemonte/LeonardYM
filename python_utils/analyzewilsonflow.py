#!/usr/bin/python

import re
import os
import sys
from subprocess import call

if len(sys.argv) != 3:
	print "This script must be called with three arguments!"
	print "Usage:\n   parallelstarter.py scriptfile folder"
	sys.exit(1)

index = 1

for namefile in os.listdir(sys.argv[2]):
	pattern = re.search('0:([0-9]+):.+',namefile)
	if pattern != None:
		tmp = int(pattern.group(1))
		call(["cp", sys.argv[1], sys.argv[1]+str(index).zfill(2)+".cfg"])
		configsfile = open(sys.argv[1]+str(index).zfill(2)+".cfg",'r')
		configsdata = configsfile.read()
		replacerule = re.search("(read_start_number\ *=\ *[0-9]+)",configsdata)
		configsdata = configsdata.replace(replacerule.group(1),"read_start_number = "+pattern.group(1))
		
		basename = re.search('output_name\ *=\ *(.+_)[0-9]+',configsdata).group(1)
		replacerule = re.search("(output_name\ *=\ *.+_[0-9]+)",configsdata)
		configsdata = configsdata.replace(replacerule.group(1),"output_name = "+basename+pattern.group(1))
		
		inputfolder = re.search('(input_directory_configurations\ *=\ *.+)',configsdata).group(1)
		configsdata = configsdata.replace(inputfolder,"input_directory_configurations = "+sys.argv[2])
		configsfile.close()
		
		outputfile = open(sys.argv[1]+str(index).zfill(2)+".cfg",'w')
		outputfile.write(configsdata)
		outputfile.close()	
		
		index = index + 1

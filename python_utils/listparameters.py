#!/usr/bin/python

import re
import os
import sys
from subprocess import call

if len(sys.argv) != 4:
	print "This script must be called with three arguments!"
	print "Usage:\n   listparameters.py folder Ns Nt"
	sys.exit(1)
	
parameters = []

for namefile in os.listdir(sys.argv[1]):
	pattern = re.search('.+\_'+sys.argv[2]+'c'+sys.argv[3]+'\_([0-9]+)b\_([0-9]+)k',namefile)
	if (pattern != None):
		if ([pattern.group(1),pattern.group(2)]) not in parameters:
			parameters.append([pattern.group(1),pattern.group(2)])

f = open(sys.argv[1]+"znaj.txt", 'w')
f.write(str(parameters).replace("[","{").replace("]","}").replace("'",""))
f.close()
	

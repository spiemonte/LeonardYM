#!/usr/bin/python

import re
import os
import sys
from subprocess import call

if len(sys.argv) != 3:
	print "This script must be called with three arguments!"
	print "Usage:\n   removeconfigurations.py folder skip"
	sys.exit(1)
	
skip = int(sys.argv[2])

for namefile in os.listdir(sys.argv[1]):
	pattern = re.search('0:([0-9]+):.+',namefile)
	if pattern != None:
		tmp = int(pattern.group(1))
		if (tmp % skip) != 0:
			call(["rm", sys.argv[1]+namefile])

for namefile in os.listdir(sys.argv[1]):
	pattern = re.search('.+\_([0-9]+)\_[0-9]+\.txt',namefile)
	if pattern != None:
		tmp = int(pattern.group(1))
		if (tmp % skip) != 0:
			call(["rm", sys.argv[1]+namefile])
	pattern = re.search('.+\_([0-9]+)\.descriptor\.txt',namefile)
	if pattern != None:
		tmp = int(pattern.group(1))
		if (tmp % skip) != 0:
			call(["rm", sys.argv[1]+namefile])

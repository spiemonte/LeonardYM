#!/usr/bin/python

import re
import os
import sys
from subprocess import call

if len(sys.argv) != 3:
	print "This script must be called with three arguments!"
	print "Usage:\n   parallelstarter.py scriptfile parallelconfigs"
	sys.exit(1)
	
f = open(sys.argv[2], 'r')
pconfigs = f.readlines()
f.close()

fields = re.split(' |\n|\t|\r',pconfigs[0])
for i in range(0,fields.count('')):
	fields.remove('')


def betaChanger(filename,beta):
	call(["./changebeta.py", filename, beta])

def kappaChanger(filename,kappa):
	call(["./changekappa.py", filename, kappa])

def NsChanger(filename,Ns):
	call(["./changeNs.py", filename, Ns])

def NtChanger(filename,Nt):
	call(["./changeNt.py", filename, Nt])
	
options = {	'Ns' : NsChanger,
			'Nt' : NtChanger,
			'beta' : betaChanger,
			'kappa' : kappaChanger
}

print "The following items will be iterated:"
for field in fields:
	print "\t"+field

pconfigs.remove(pconfigs[0])

index = 1

for line in pconfigs:
	if re.search('[0-9]',line) != None:
		call(["cp", sys.argv[1], sys.argv[1]+str(index).zfill(2)+".cfg"])
		values = re.split(' |\n|\t|\r',line)
		for i in range(0,values.count('')):
			values.remove('')
		if (len(values) != len(fields)):
			print "Line \""+line.replace('\n','')+"\" off of range!"
		for field, value in zip(fields,values):
			options[field](sys.argv[1]+str(index).zfill(2)+".cfg",value)
		index = index + 1
	



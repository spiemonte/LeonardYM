#!/usr/bin/python
import os
import sys
import re
import numpy
import getopt
import random
import math
import time
import subprocess

if len(sys.argv) != 3:
	print "This script must be called with two argument!"
	print "Usage:\n   simulation_generator.py generatorfile configfile"
	sys.exit(1)
	
inputfile = open(sys.argv[1],'r')
filedata = inputfile.read()

def get_option(configsname):
	confspattern = re.search("^"+configsname+"\ *=\ *([^\n]+)", filedata, re.MULTILINE)
	if confspattern != None:
		result = re.split(",\ *", confspattern.group(1))
		return result
	else:
		return None

number_simulations = int(get_option("number_simulations")[0])

for i in range(1,number_simulations+1):
	subprocess.call("cp "+sys.argv[2]+" "+sys.argv[2]+("{0:02d}".format(i))+".cfg",shell=True)
	
def set_option(optname,commandname):
	opt = get_option(optname)
	if opt != None:
		if len(opt) != number_simulations:
			if len(opt) == 1:
				for i in range(1,number_simulations+1):
					print commandname+sys.argv[2]+("{0:02d}".format(i))+".cfg "+opt[0]
					subprocess.call(commandname+sys.argv[2]+("{0:02d}".format(i))+".cfg "+opt[0],shell=True)
			else:
				print "wrong number of ", optname, "(", len(opt),") provided"
		else:
			for pair in zip(opt,range(1,number_simulations+1)):
				print commandname+sys.argv[2]+("{0:02d}".format(pair[1]))+".cfg "+pair[0]
				subprocess.call(commandname+sys.argv[2]+("{0:02d}".format(pair[1]))+".cfg "+pair[0],shell=True)
			
set_option("nt", "./changeNt.py ")
set_option("ns", "./changeNs.py ")
set_option("beta", "./changebeta.py ")
set_option("kappa", "./changekappa.py ")


#!/usr/bin/python
import os
import sys
import re
import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import getopt
import random
import math

folder = ""
basename = ""
init=100
blocking=3

def usage():
	print "Usage: hmc_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: hmc_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking]"

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:b:s:k:", ["help", "folder=","basename=","skip=","blocking="])
except getopt.GetoptError, err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)
for option, argument in opts:
	if option in ("-h", "--help"):
		usage()
	elif option in ("-f", "--folder"):
		folder = argument
	elif option in ("-b", "--basename"):
		basename = argument
	elif option in ("-s", "--skip"):
		init = int(argument)
	elif option in ("-k", "--blocking"):
		blocking = int(argument)
if basename == "":
	print "Basename not set!"
	usage()
	sys.exit(2)
if folder == "":
	print "Folder not set!"
	usage()
	sys.exit(2)
	
print "Skipping "+str(init)+" initial configurations as thermalization"
print "Using "+str(blocking)+" blocking levels"

import load_vector

data = load_vector.load(basename,"hmc_history.txt",folder)
finaldata = [j for i in data for j in i]

end=len(finaldata)

if end == 0:
	print "Data are not correctly loaded, do they exist?"
	sys.exit(2)

toplot = zip(*finaldata)

import bootstrapping as boot

title = "HMC history for run "+basename

acceptance = numpy.mean(toplot[1][init:end])
resultstr = "Acceptance rate: "+str(round(acceptance,4))
print resultstr
title += "\n"+resultstr

def erf(x):
	# save the sign of x
	sign = 1 if x >= 0 else -1
	x = abs(x)

	# constants
	a1 =  0.254829592
	a2 = -0.284496736
	a3 =  1.421413741
	a4 = -1.453152027
	a5 =  1.061405429
	p  =  0.3275911

	# A&S formula 7.1.26
	t = 1.0/(1.0 + p*x)
	y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
	return sign*y # erf(-x) = -erf(x)

creutz_equality = boot.bootstrap(boot.block([math.exp(-i) for i in toplot[0][init:end]],blocking))
roundlevel = int(-math.log10(creutz_equality[1]) + 2)
resultstr = "Creutz equality: "+str(round(creutz_equality[0],roundlevel))+" +/- "+str(round(creutz_equality[1],roundlevel))
print resultstr
title += "\n"+resultstr

expected_acceptance = 1. - erf(math.sqrt(numpy.mean(toplot[0][init:end])/2.))
resultstr = "Expected acceptance: "+str(round(expected_acceptance,4))
print resultstr
title += "\n"+resultstr

plt.suptitle(title)

ax = plt.subplot(2,1,1)
ax.plot(range(init,end),toplot[0][init:end],color='b')
ax.yaxis.set_label_text("Energy variation")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.xaxis.set_major_locator(plticker.LinearLocator(5))

ax.xaxis.set_label_text("T Monte Carlo")

ax = plt.subplot(2,1,2)
ax.hist(toplot[0][init:end],facecolor='b')
ax.xaxis.set_major_locator(plticker.LinearLocator(5))
ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.yaxis.set_label_text("Frequency")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

ax.xaxis.set_label_text("Energy variation value")

plt.subplots_adjust(wspace=0.3,left=0.2,top=0.83,bottom=0.1)
plt.savefig(basename+"_hmc_history")





#!/usr/bin/python
import os
import sys
import re
import numpy
import getopt
import random
import math

folder = ""
basename = ""
init=100
blocking=3

def usage():
	print "Usage: chiral_condensate_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: chiral_condensate_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking]"

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

data = load_vector.load(basename,"condensate.txt",folder)
datac = load_vector.load(basename,"condensate_connected.txt",folder)

finaldata = [j for i in data for j in i]
finaldatac = [j for i in datac for j in i]

finaldata_cleaned = []
finaldatac_cleaned = []
removed = 0

for i in range(len(finaldata)):
	if finaldata[i][1] > 0.1:
		removed = removed + 1
	else:
		finaldata_cleaned.append(finaldata[i])
		finaldatac_cleaned.append(finaldatac[i])
		
finaldata = finaldata_cleaned
finaldatac = finaldatac_cleaned

end=len(finaldata)

toplot = zip(*finaldata)

print "Removing ", removed, " exceptional configurations."

import bootstrapping as boot

title = "Chiral condensate for run "+basename

chiral_condensate = boot.bootstrap(boot.block(toplot[0][init:end],blocking))

chiral_condensate_susc = boot.bootstrap_susc(boot.block([abs(toplot[0][i]) for i in range(init,end)],blocking),boot.block([toplot[0][i]*toplot[0][i] for i in range(init,end)],blocking))

connected = zip(*finaldatac)
chiral_condensate_connected_susc = boot.bootstrap(boot.block(connected[0][init:end],blocking))


roundlevel = int(-math.log10(chiral_condensate[1]) + 2)
resultstr = "Chiral condensate expectation value: "+str(round(chiral_condensate[0],roundlevel))+" +/- "+str(round(chiral_condensate[1],roundlevel))

roundlevel_susc = int(-math.log10(chiral_condensate_susc[1]) + 2)
resultstr += "\nChiral condensate disconnected susceptibility: "+str(round(chiral_condensate_susc[0],roundlevel_susc))+" +/- "+str(round(chiral_condensate_susc[1],roundlevel_susc))

roundlevel_connected_susc = int(-math.log10(chiral_condensate_connected_susc[1]) + 2)
resultstr += "\nChiral condensate connected susceptibility: "+str(round(chiral_condensate_connected_susc[0],roundlevel_connected_susc))+" +/- "+str(round(chiral_condensate_connected_susc[1],roundlevel_connected_susc))

print resultstr
title += "\n"+resultstr

pattern = re.search('([0-9]+)c([0-9]+)',basename)
latticedata = ""
if pattern != None:
	latticedata += pattern.group(1)+", "+pattern.group(2)+", "
pattern = re.search('([0-9]+)b',basename)
if pattern != None:
	latticedata += pattern.group(1)[0]+"."+pattern.group(1)[1:]+", "
pattern = re.search('([0-9]+)k',basename)
if pattern != None:
	latticedata += "0."+pattern.group(1)+", "

print "Mathematica output: {"+latticedata+(str(round(chiral_condensate[0],roundlevel))+", "+str(round(chiral_condensate[1],roundlevel))+", "+str(round(chiral_condensate_susc[0],roundlevel_susc))+", "+str(round(chiral_condensate_susc[1],roundlevel_susc))+", "+str(round(chiral_condensate_connected_susc[0],roundlevel_connected_susc))+", "+str(round(chiral_condensate_connected_susc[1],roundlevel_connected_susc))).replace("e","*10^")+"}"

try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.ticker as plticker

	plt.suptitle(title)

	ax = plt.subplot(1,2,1)
	ax.plot(range(init,end),toplot[0][init:end],color='b')
	ax.yaxis.set_label_text("Chiral\ncondensate")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))

	ax.xaxis.set_label_text("T Monte Carlo")

	ax = plt.subplot(1,2,2)
	ax.hist(toplot[0][init:end],facecolor='b')
	ax.xaxis.set_major_locator(plticker.LinearLocator(3))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
	ax.yaxis.set_label_text("Frequency")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

	ax.xaxis.set_label_text("Chiral condensate value")

	plt.subplots_adjust(wspace=0.3,left=0.2,top=0.83,bottom=0.1)
	plt.savefig(basename+"_chiral_condensate")
except ImportError:
	print "Matplotlib not present, no plot output will be produced!"





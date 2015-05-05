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
	print "Usage: plaquette_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: plaquette_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking]"

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

data = load_vector.load(basename,"plaquette.txt",folder)
finaldata = [j for i in data for j in i]

end=len(finaldata)

if end == 0:
	print "Data are not correctly loaded, do they exist?"
	sys.exit(2)
elif end < init:
	print "Not enough data to skip the first ",init," configs!"
	print "Use the option -s to a lower value than ", end
	sys.exit(2)
else:
	print "Analysing ", end-init, " configurations"

toplot = zip(*finaldata)

import bootstrapping as boot
import utils

title = "Plaquette for run "+basename

total_plaquette = boot.bootstrap(boot.block(toplot[0][init:end],blocking))
resultstr = "Total plaquette expectation value: "+utils.str_round(total_plaquette[0],total_plaquette[1])+" +/- "+utils.str_round(total_plaquette[1],total_plaquette[1])
print resultstr
title += "\n"+resultstr

spatial_plaquette = boot.bootstrap(boot.block(toplot[1][init:end],blocking))
resultstr = "Spatial plaquette expectation value: "+utils.str_round(spatial_plaquette[0],spatial_plaquette[1])+" +/- "+utils.str_round(spatial_plaquette[1],spatial_plaquette[1])
print resultstr
title += "\n"+resultstr

temporal_plaquette = boot.bootstrap(boot.block(toplot[2][init:end],blocking))
resultstr = "Temporal plaquette expectation value: "+utils.str_round(temporal_plaquette[0],temporal_plaquette[1])+" +/- "+utils.str_round(temporal_plaquette[1],temporal_plaquette[1])
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
	
print "Mathematica output: {"+latticedata+(utils.str_round(total_plaquette[0],total_plaquette[1])+", "+utils.str_round(total_plaquette[1],total_plaquette[1])+", "+utils.str_round(spatial_plaquette[0],spatial_plaquette[1])+", "+utils.str_round(spatial_plaquette[1],spatial_plaquette[1])+", "+utils.str_round(temporal_plaquette[0],temporal_plaquette[1])+", "+utils.str_round(temporal_plaquette[1],temporal_plaquette[1])).replace("e","*10^")+"}"

try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.ticker as plticker

	plt.suptitle(title)

	ax = plt.subplot(3,2,1)
	ax.plot(range(init,end),toplot[0][init:end],color='b')
	ax.yaxis.set_label_text("Total\nplaquette")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))

	ax = plt.subplot(3,2,3)
	ax.plot(range(init,end),toplot[1][init:end],color='r')
	ax.yaxis.set_label_text("Spatial\nplaquette")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))

	ax = plt.subplot(3,2,5)
	ax.plot(range(init,end),toplot[2][init:end],color='g')
	ax.yaxis.set_label_text("Temporal\nplaquette")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))

	ax.xaxis.set_label_text("T Monte Carlo")

	ax = plt.subplot(3,2,2)
	ax.hist(toplot[0][init:end],facecolor='b')
	ax.xaxis.set_major_locator(plticker.LinearLocator(3))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
	ax.yaxis.set_label_text("Frequency")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

	ax = plt.subplot(3,2,4)
	ax.hist(toplot[1][init:end],facecolor='r')
	ax.xaxis.set_major_locator(plticker.LinearLocator(3))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
	ax.yaxis.set_label_text("Frequency")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

	ax = plt.subplot(3,2,6)
	ax.hist(toplot[2][init:end],facecolor='g')
	ax.xaxis.set_major_locator(plticker.LinearLocator(3))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
	ax.yaxis.set_label_text("Frequency")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

	ax.xaxis.set_label_text("Plaquette value")

	plt.subplots_adjust(wspace=0.3,left=0.2,top=0.83,bottom=0.1)
	plt.savefig(basename+"_plaquette")
except ImportError:
	print "Matplotlib not present, no plot output will be produced!"






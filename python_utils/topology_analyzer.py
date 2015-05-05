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
number_boot_resamplings = 100
reference = 0.3
tr = 0.
dd = 1

def usage():
	print "Usage: wilson_flow_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: wilson_flow_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking --number_boot_resamplings=n --top_reference=tr --scale_reference=r]"

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:b:s:k:", ["help", "folder=","basename=","skip=","blocking=","number_boot_resamplings=","top_reference=","scale_reference=","delta_derivative="])
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
	elif option == "--number_boot_resamplings":
		number_boot_resamplings = int(argument)
	elif option == "--scale_reference":
		reference = float(argument)
	elif option == "--top_reference":
		tr = float(argument)
if basename == "":
	print "Basename not set!"
	usage()
	sys.exit(2)
if folder == "":
	print "Folder not set!"
	usage()
	sys.exit(2)
	
print "Skipping ", init, " initial configurations as thermalization"
print "Using ", blocking," blocking levels"
print "The scale will be set at the reference value ", reference, " (use --scale_reference to change it)"
print "The topological charge will be measured at the flow time ", tr, " (use --top_reference to change it)"
print "The error will be estimated with ", number_boot_resamplings, " bootstrapping resamplings (use --number_boot_resamplings to change it)"

import load_vector

data_wilson_flow = load_vector.load(basename,"wilson_flow.txt",folder)
finaldata_wilson_flow = [j for i in data_wilson_flow for j in i]

data_topological_charge = load_vector.load(basename,"topological_charge.txt",folder)
finaldata_topological_charge = [j for i in data_topological_charge for j in i]

if len(finaldata_wilson_flow) != len(finaldata_topological_charge):
	print "Warning, mismatch between the size of the data of the topological charge and of the scale!"
	finaldata_topological_charge = finaldata_topological_charge[:min(len(finaldata_wilson_flow),len(finaldata_topological_charge))]
	finaldata_wilson_flow = finaldata_wilson_flow[:min(len(finaldata_wilson_flow),len(finaldata_topological_charge))]

end=len(finaldata_wilson_flow)

if end == 0:
	print "Data are not correctly loaded, do they exist?"
	sys.exit(2)
elif end < init:
	print "Not enough data to skip the first ",init," configs!"
	print "Use the option -s to a lower value than ", end
	sys.exit(2)
else:
	print "Analysing ", end-init, " configurations"

toplotwf = zip(*finaldata_wilson_flow)
toplotq = zip(*finaldata_topological_charge)

print "The flow has", len(toplotwf), "discretization intervals, to the final flow time tau =", toplotwf[-1][0][0]

if toplotwf[-1][0][0] < tr:
	print "Time ", tr, " to measure the topological charge larger than ", toplotwf[-1][0][0], " too big! (use --top_reference to change it)"
	sys.exit(2)



import bootstrapping as boot
import optimize as opt

wilson_flow_data = [ [ttenergy[0][0]] + boot.bootstrap(boot.block(zip(*ttenergy)[1][init:end],blocking),50) for ttenergy in toplotwf]
topological_susceptibility_data = [ [ttq[0][0]] + boot.bootstrap_susc(boot.block(zip(*ttq)[1][init:end],blocking), boot.block([ttq[i][1]*ttq[i][1] for i in range(init,end)],blocking),50) for ttq in toplotq]


def auto_w0(wilson_flow_data_boot):
	for i in range(dd,len(wilson_flow_data_boot)-dd):
		derivative2 = wilson_flow_data_boot[i][0]*(wilson_flow_data_boot[i+dd][1]-wilson_flow_data_boot[i-dd][1])/(wilson_flow_data_boot[i+dd][0]-wilson_flow_data_boot[i-dd][0])
		if derivative2 > reference:
			derivative1 = wilson_flow_data_boot[i-1][0]*(wilson_flow_data_boot[i+dd-1][1]-wilson_flow_data_boot[i-dd-1][1])/(wilson_flow_data_boot[i+dd-1][0]-wilson_flow_data_boot[i-dd-1][0])
			t0 = wilson_flow_data_boot[i-1][0]
			delta = wilson_flow_data_boot[i][0] - wilson_flow_data_boot[i-1][0]
			return math.sqrt(t0 + delta*(reference - derivative1)/(derivative2 - derivative1))
	print "Warning, auto w0 calculation out of range!"
	return math.sqrt(wilson_flow_data_boot[-1][0])
			
def auto_t0(wilson_flow_data_boot):
	for i in range(1,len(wilson_flow_data_boot)):
		if wilson_flow_data_boot[i][1] > reference:
			t0 = wilson_flow_data_boot[i-1][0]
			delta = wilson_flow_data_boot[i][0] - wilson_flow_data_boot[i-1][0]
			return t0 + delta*(reference - wilson_flow_data_boot[i-1][1])/(wilson_flow_data_boot[i][1] - wilson_flow_data_boot[i-1][1])
	print "Warning, auto t0 calculation out of range!"
	return wilson_flow_data_boot[-1][0]
	
def topological_charge(qdata):
	for value in zip(qdata, qdata[1:]):
		if value[1][0] > tr:
			return ((value[1][1]-value[0][1])/(value[1][0]-value[0][0]))*(tr-value[0][0])+value[0][1]
	
def measure_susceptibility(susceptibility_data):
	for value in zip(susceptibility_data, susceptibility_data[1:]):
		if value[1][0] > tr:
			return [((value[1][1]-value[0][1])/(value[1][0]-value[0][0]))*(tr-value[0][0])+value[0][1],((value[1][2]-value[0][2])/(value[1][0]-value[0][0]))*(tr-value[0][0])+value[0][2]]
	
print "Wilson flow scale t_", reference, ":" , auto_t0(wilson_flow_data)
susceptibility = measure_susceptibility(topological_susceptibility_data)

print "Topological susceptibility at t=", tr, ":" ,susceptibility[0], " +/- ", susceptibility[1]

print "Running bootstrapping analysis ..."
	
def bootstrap_topological_susceptibility(datawf,dataq,datasq,trials=500):
	boot_data = []
	orig = [(auto_t0(wilson_flow_data)**2)*measure_susceptibility(topological_susceptibility_data)[0],auto_t0(wilson_flow_data)]
	for step in range(1,trials):
		randomset = [
			random.randint(0,len(datawf[0][1])-1)
				for i in range(len(datawf[0][1]))
			]
		wilson_flow_data_boot = [[ttenergy[0], numpy.mean([ttenergy[1][j] for j in randomset]), 0.] for ttenergy in datawf]
		topological_susceptibility_data_boot = [ [dataq[i][0], numpy.mean([datasq[i][1][j] for j in randomset]) - (numpy.mean([dataq[i][1][j] for j in randomset])**2) , 0.] for i in range(len(dataq))]
		
		t0 = auto_t0(wilson_flow_data_boot)
		boot_data.append([(t0**2)*measure_susceptibility(topological_susceptibility_data_boot)[0],t0])
	return [[2*orig[i] - numpy.mean(zip(*boot_data)[i]),math.sqrt(float(len(datawf[0])-1)/len(datawf[0]))*numpy.std(zip(*boot_data)[i])] for i in range(2)]
	

physical_topological_susceptibility = bootstrap_topological_susceptibility(
	[[ttenergy[0][0], boot.block(zip(*ttenergy)[1][init:end],blocking)] for ttenergy in toplotwf],
	[[ttq[0][0], boot.block(zip(*ttq)[1][init:end],blocking)] for ttq in toplotq],
	[[ttq[0][0], boot.block([ttq[i][1]*ttq[i][1] for i in range(init,end)],blocking)] for ttq in toplotq],
	number_boot_resamplings)

title = "Topological susceptibility for run "+basename

import utils

resultstr = "Topological susceptibility chi: "+utils.str_round(susceptibility[0],susceptibility[1])+" +/- "+utils.str_round(susceptibility[1],susceptibility[1])
print resultstr
title += "\n"+resultstr

resultstr = "Scale t_"+str(reference)+": "+utils.str_round(physical_topological_susceptibility[1][0],physical_topological_susceptibility[1][1])+" +/- "+utils.str_round(physical_topological_susceptibility[1][1],physical_topological_susceptibility[1][1])
print resultstr
title += "\n"+resultstr

resultstr = "Topological susceptibility chi*(t_"+str(reference)+")^2: "+utils.str_round(physical_topological_susceptibility[0][0],physical_topological_susceptibility[0][1])+" +/- "+utils.str_round(physical_topological_susceptibility[0][1],physical_topological_susceptibility[0][1])
print resultstr
title += "\n"+resultstr

pattern = re.search('([0-9]+)c([0-9]+)',basename)
latticedata = ""
moutput = ""
if pattern != None:
	latticedata += pattern.group(1)+", "+pattern.group(2)+", "
	volume = (float(pattern.group(1))**3)*float(pattern.group(2))
	resultstr = "Topological susceptibility chi*(t_"+str(reference)+")^2/V: "+utils.str_round(physical_topological_susceptibility[0][0]/volume,physical_topological_susceptibility[0][1]/volume)+" +/- "+utils.str_round(physical_topological_susceptibility[0][1]/volume,physical_topological_susceptibility[0][1]/volume)
	print resultstr
	title += "\n"+resultstr
	
	moutput = (
		utils.str_round(physical_topological_susceptibility[0][0]/volume,physical_topological_susceptibility[0][1]/volume)
		+", "+utils.str_round(physical_topological_susceptibility[0][1]/volume,physical_topological_susceptibility[0][1]/volume)
		+", "+utils.str_round(physical_topological_susceptibility[1][0],physical_topological_susceptibility[1][1])
		+", "+utils.str_round(physical_topological_susceptibility[1][1],physical_topological_susceptibility[1][1])).replace("e","*10^")
else:
	moutput = (
		utils.str_round(physical_topological_susceptibility[0][0],physical_topological_susceptibility[0][1])
		+", "+utils.str_round(physical_topological_susceptibility[0][1],physical_topological_susceptibility[0][1])
		+", "+utils.str_round(physical_topological_susceptibility[1][0],physical_topological_susceptibility[1][1])
		+", "+utils.str_round(physical_topological_susceptibility[1][1],physical_topological_susceptibility[1][1])).replace("e","*10^")
	
pattern = re.search('([0-9]+)b',basename)
if pattern != None:
	latticedata += pattern.group(1)[0]+"."+pattern.group(1)[1:]+", "
pattern = re.search('([0-9]+)k',basename)
if pattern != None:
	latticedata += "0."+pattern.group(1)+", "
	
print "Mathematica output: {"+latticedata+moutput+"}"

try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.ticker as plticker

	plt.suptitle(title)

	ax = plt.subplot(3,1,1)
	ax.errorbar([elem[0] for elem in topological_susceptibility_data],[elem[1] for elem in topological_susceptibility_data],color='b',yerr=[elem[2] for elem in wilson_flow_data], marker='o',linestyle='')
	ax.yaxis.set_label_text("chi")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.2f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, toplotwf[-1][0][0])
	
	ax.xaxis.set_label_text("t")
	
	history = [topological_charge(topological_charge_flow) for topological_charge_flow in finaldata_topological_charge]
			
	ax = plt.subplot(3,1,2)
	ax.plot(range(end-init),history[init:end],color='b')
	ax.yaxis.set_label_text("Q_top")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.2f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, end-init)

	ax.xaxis.set_label_text("Configuration")
	
	ax = plt.subplot(3,1,3)
	ax.hist(history[init:end],facecolor='b',bins=numpy.arange(min(history[init:end]), max(history[init:end]) + 0.1, 0.1))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.1f'))
	ax.yaxis.set_label_text("Frequency")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

	ax.xaxis.set_label_text("Q_top")

	plt.subplots_adjust(hspace=0.45,wspace=0.35,left=0.2,top=0.81,bottom=0.1)
	plt.savefig(basename+"_topological_charge")

except ImportError:
	print "Matplotlib not present, no plot output will be produced!"

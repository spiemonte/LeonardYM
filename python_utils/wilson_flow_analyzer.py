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
measure_coupling = False
tc = 0.
delta= 3

order=1

def usage():
	print "Usage: wilson_flow_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: wilson_flow_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking --number_boot_resamplings=n --reference=r --order=n --t_coupling=tc --delta=d]"

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:b:s:k:", ["help", "folder=","basename=","skip=","blocking=","number_boot_resamplings=","t_coupling=","reference=","order=","delta="])
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
	elif option == "--reference":
		reference = float(argument)
	elif option == "--order":
		order = int(argument)
	elif option == "--delta":
		delta = int(argument)
	elif option == "--t_coupling":
		measure_coupling = True
		tc = float(argument)
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
print "The scale will be set at the reference value ", reference, " (use --reference to change it)"
print "A polynomial lest square of order ", order, " will be used (use --order to change it)"
print " fitting ", 2*delta + 1, " points (use --delta to change it)"
print "The error will be estimated with ", number_boot_resamplings, " bootstrapping resamplings (use --number_boot_resamplings to change it)"

import load_vector

data_wilson_flow = load_vector.load(basename,"wilson_flow.txt",folder)
finaldata_wilson_flow = [j for i in data_wilson_flow for j in i]

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

print "The flow has", len(toplotwf), "discretization intervals, to the final flow time tau =", toplotwf[-1][0][0]
print "The discretization interval of the flow is ", toplotwf[-1][0][0]/len(toplotwf)

if measure_coupling:
	if tc > toplotwf[-1][0][0]:
		print "Time", tc, "used to specify the coupling too large (use --t_coupling to change it)"
		measure_coupling = False

import bootstrapping as boot
import optimize as opt

wilson_flow_data = [ [ttenergy[0][0], numpy.mean(zip(*ttenergy)[1][init:end])] for ttenergy in toplotwf]

#wfdata = [[(wilson_flow_data[t][0],),wilson_flow_data[t][1],wilson_flow_data[t][2]] for t in range(td,tf)]

#trial, minimum = opt.parallel_tempering(chibtWf, [0.3,5.])

#resultwf = opt.minimum_chisq(thWf,wfdata,tuple(trial),derivative_function=d_thWf)

#print "Wilson flow scale t_", reference, ":" , (reference - resultwf[1])/resultwf[0]
#print "Wilson flow scale w_", reference, ":" , math.sqrt((reference)/resultwf[0])

print "Running bootstrapping analysis ..."
#chisq = chibtWf(*resultwf)/(tf-td-3)

import utils

def auto_w0(wilson_flow_data_boot):
	result = wilson_flow_data_boot[-1][0]
	for i in range(1+delta,len(wilson_flow_data_boot)-delta):
		points = [wilson_flow_data_boot[j][0] for j in range(i-delta,i+delta+1)]
		values = [wilson_flow_data_boot[j][1] for j in range(i-delta,i+delta+1)]
		derivative = utils.Polynomial([0,1])*utils.LeastSquarePolynomial(points,values,order).derive()+utils.Polynomial([-reference])
		if derivative.evaluate(points[0])*derivative.evaluate(points[-1]) < 0:
			result = derivative.zero(points[0],points[-1]);
			if result > points[len(points)/2-1]:
				return math.sqrt(result)
	print "Warning, auto w0 calculation out of range!"
	return math.sqrt(result)

def auto_t0(wilson_flow_data_boot):
	result = wilson_flow_data_boot[-1][0]
	for i in range(1+delta,len(wilson_flow_data_boot)-delta):
		points = [wilson_flow_data_boot[j][0] for j in range(i-delta,i+delta+1)]
		values = [wilson_flow_data_boot[j][1] for j in range(i-delta,i+delta+1)]
		function = utils.LeastSquarePolynomial(points,values,order)+utils.Polynomial([-reference])
		if function.evaluate(points[0])*function.evaluate(points[-1]) < 0:
			result = function.zero(points[0],points[-1]);
			if result > points[len(points)/2-1]:
				return result
	print "Warning, auto t0 calculation out of range!"
	return result
			
#def auto_t0(wilson_flow_data_boot):
#	for i in range(1,len(wilson_flow_data_boot)):
#		if wilson_flow_data_boot[i][1] > reference:
#			t0 = wilson_flow_data_boot[i-1][0]
#			delta = wilson_flow_data_boot[i][0] - wilson_flow_data_boot[i-1][0]
#			return t0 + delta*(reference - wilson_flow_data_boot[i-1][1])/(wilson_flow_data_boot[i][1] - wilson_flow_data_boot[i-1][1])
#	print "Warning, auto t0 calculation out of range!"
#	return wilson_flow_data_boot[-1][0]
	
	
	
def bootstrap_scale(datawf,trials=500):
	boot_data = []
	
	orig = [auto_t0(wilson_flow_data), auto_w0(wilson_flow_data)]
	for step in range(1,trials):
		randomset = [random.randint(0,len(datawf[0][1])-1) for i in range(len(datawf[0][1]))]
		wilson_flow_data_boot = [[ttenergy[0], numpy.mean([ttenergy[1][j] for j in randomset])] for ttenergy in datawf]		
		
		boot_data.append([auto_t0(wilson_flow_data_boot), auto_w0(wilson_flow_data_boot)])
	return [[2*orig[i] - numpy.mean(zip(*boot_data)[i]),math.sqrt(float(len(datawf[0])-1)/len(datawf[0]))*numpy.std(zip(*boot_data)[i])] for i in range(2)]
	

wilson_fit = bootstrap_scale([[ttenergy[0][0], boot.block(zip(*ttenergy)[1][init:end],blocking)] for ttenergy in toplotwf],number_boot_resamplings)

title = "Wilson flow scale for run "+basename

resultstr = "Auto scale t_"+str(reference)+": "+utils.str_round(wilson_fit[0][0],wilson_fit[0][1])+" +/- "+utils.str_round(wilson_fit[0][1],wilson_fit[0][1])
print resultstr
title += "\n"+resultstr

resultstr = "Auto scale w_"+str(reference)+": "+utils.str_round(wilson_fit[1][0],wilson_fit[1][1])+" +/- "+utils.str_round(wilson_fit[1][1],wilson_fit[1][1])
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

def coupling(wilson_flow_data_boot,tt):
        result = wilson_flow_data_boot[-1][1]
        for i in range(1+delta,len(wilson_flow_data_boot)-delta):
		if wilson_flow_data_boot[i][0] > tc:
                	points = [wilson_flow_data_boot[j][0] for j in range(i-delta,i+delta+1)]
                	values = [wilson_flow_data_boot[j][1] for j in range(i-delta,i+delta+1)]
                	function = utils.LeastSquarePolynomial(points,values,order)
			return function.evaluate(tt)
        print "Warning, coupling calculation out of range!"
        return result

def bootstrap_coupling(datawf,trials=500):
	boot_data = []

        orig = coupling(wilson_flow_data,tc)
        for step in range(1,trials):
                randomset = [random.randint(0,len(datawf[0][1])-1) for i in range(len(datawf[0][1]))]
                wilson_flow_data_boot = [[ttenergy[0], numpy.mean([ttenergy[1][j] for j in randomset])] for ttenergy in datawf]

                boot_data.append(coupling(wilson_flow_data_boot,tc))
        return [2*orig - numpy.mean(boot_data),math.sqrt(float(len(datawf[0])-1)/len(datawf[0]))*numpy.std(boot_data)]
	
coupling_str = ""
if measure_coupling:
	cc = bootstrap_coupling([[ttenergy[0][0], boot.block(zip(*ttenergy)[1][init:end],blocking)] for ttenergy in toplotwf],number_boot_resamplings)

	print "Gauge coupling defined at ", tc, ": ", utils.str_round(cc[0],cc[1]), " +/- ", utils.str_round(cc[1],cc[1])
	coupling_str = ", "+utils.str_round(cc[0],cc[1])+", "+utils.str_round(cc[1],cc[1])
	
			

print "Mathematica output: {"+latticedata+(utils.str_round(wilson_fit[0][0],wilson_fit[0][1])+", "+utils.str_round(wilson_fit[0][1],wilson_fit[0][1])+", "+utils.str_round(wilson_fit[1][0],wilson_fit[1][1])+", "+utils.str_round(wilson_fit[1][1],wilson_fit[1][1])+coupling_str).replace("e","*10^")+"}"

try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.ticker as plticker

	plt.suptitle(title)

	ax = plt.subplot(2,1,1)
	ax.plot([elem[0] for elem in wilson_flow_data],[elem[1] for elem in wilson_flow_data],color='b')
	ax.plot(numpy.arange(0.,wilson_fit[0][0],0.02),[reference for elem in numpy.arange(0.,wilson_fit[0][0],0.02)],color='g')
	ax.plot([wilson_fit[0][0] for elem in numpy.arange(0.,reference,0.01)],numpy.arange(0.,reference,0.01),color='g')
	#ax.plot([elem[0] for elem in wilson_flow_data],[thWf(elem[0],wilson_fit[0][0],wilson_fit[1][0]) for elem in wilson_flow_data],color='r')
	ax.yaxis.set_label_text("t*t*<E>")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.2f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, toplotwf[-1][0][0])
	
	ax.xaxis.set_label_text("t")
	
	history = []
	for config in finaldata_wilson_flow:
		t = len(history)
		for elem in config:
			if elem[1] > reference:
				history.append(math.sqrt(elem[0]))
				break
		if t == len(history):
			history.append(math.sqrt(config[-1][0]))
			
	ax = plt.subplot(2,1,2)
	ax.plot(range(end-init),history[init:end],color='b')
	ax.yaxis.set_label_text("sqrt(t_"+str(reference)+")")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.2f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, end-init)

	ax.xaxis.set_label_text("Configuration")

	plt.subplots_adjust(hspace=0.3,wspace=0.3,left=0.2,top=0.86,bottom=0.1)
	plt.savefig(basename+"_wilson_flow")

except ImportError:
	print "Matplotlib not present, no plot output will be produced!"

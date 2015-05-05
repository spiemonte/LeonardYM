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
td = 5
tf = 10
reference = 0.3
measure_coupling = False
tc = 0.

def usage():
	print "Usage: wilson_flow_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: wilson_flow_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking --number_boot_resamplings=n --first_t_fit=t --last_t_fit=t --reference=r --t_coupling=tc]"

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:b:s:k:", ["help", "folder=","basename=","skip=","blocking=","number_boot_resamplings=","first_t_fit=","last_t_fit=","t_coupling=","reference="])
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
	elif option == "--first_t_fit":
		td = int(argument)
	elif option == "--last_t_fit":
		tf = int(argument)
	elif option == "--reference":
		reference = float(argument)
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
print "The fit will start from t =", td," (use --first_t_fit to change it)"
print "The fit will end to t =", tf," (use --last_t_fit to change it)"
print "The scale will be set at the reference value ", reference, " (use --reference to change it)"
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
elif tf <= td:
	print "Initial t for the fit larger than the final!"
	sys.exit(2)
else:
	print "Analysing ", end-init, " configurations"

toplotwf = zip(*finaldata_wilson_flow)

print "The flow has", len(toplotwf), "discretization intervals, to the final flow time tau =", toplotwf[-1][0][0]
print " the fit interval is (",toplotwf[td][0][0],", ", toplotwf[tf][0][0], ")"

if measure_coupling:
	if tc > toplotwf[-1][0][0]:
		print "Time", tc, "used to specify the coupling too large (use --t_coupling to change it)"
		measure_coupling = False

import bootstrapping as boot
import optimize as opt

wilson_flow_data = [ [ttenergy[0][0]] + boot.bootstrap(boot.block(zip(*ttenergy)[1][init:end],blocking)) for ttenergy in toplotwf]

def thWf(t, m, q):
	return m*t + q

def d_thWf(t, m, q):
	return [t, 1]
	
def chibtWf(m, q):
	return sum([((thWf(wilson_flow_data[t][0], m, q) - wilson_flow_data[t][1])**2)/(wilson_flow_data[t][2]**2) for t in range(td,tf)])


wfdata = [[(wilson_flow_data[t][0],),wilson_flow_data[t][1],wilson_flow_data[t][2]] for t in range(td,tf)]

trial, minimum = opt.parallel_tempering(chibtWf, [0.3,5.])

resultwf = opt.minimum_chisq(thWf,wfdata,tuple(trial),derivative_function=d_thWf)

print "Wilson flow scale t_", reference, ":" , (reference - resultwf[1])/resultwf[0]
print "Wilson flow scale w_", reference, ":" , math.sqrt((reference)/resultwf[0])

print "Running bootstrapping analysis ..."
chisq = chibtWf(*resultwf)/(tf-td-3)

def auto_w0(wilson_flow_data_boot):
	for i in range(2,len(wilson_flow_data_boot)-1):
		derivative2 = wilson_flow_data_boot[i][0]*(wilson_flow_data_boot[i+1][1]-wilson_flow_data_boot[i-1][1])/(wilson_flow_data_boot[i+1][0]-wilson_flow_data_boot[i-1][0])
		if derivative2 > reference:
			derivative1 = wilson_flow_data_boot[i-1][0]*(wilson_flow_data_boot[i][1]-wilson_flow_data_boot[i-2][1])/(wilson_flow_data_boot[i][0]-wilson_flow_data_boot[i-2][0])
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
	
def bootstrap_scale(datawf,trials=500):
	boot_data = []
	
	origwf = opt.minimum_chisq(thWf,wfdata,resultwf,derivative_function=d_thWf)
	
	orig = [origwf[0], origwf[1], (reference - origwf[1])/origwf[0], math.sqrt((reference)/origwf[0]), auto_t0(wilson_flow_data), auto_w0(wilson_flow_data)]
	for step in range(1,trials):
		#wilson_flow_data_boot = [[ttenergy[0]]+boot.bootstrap([ttenergy[1][j] for j in [random.randint(0,len(ttenergy[1])-1) for i in range(len(ttenergy[1]))]]) for ttenergy in datawf]
		#wilson_flow_data_boot = [[ttenergy[0], numpy.mean([ttenergy[1][j] for j in [random.randint(0,len(ttenergy[1])-1) for i in range(len(ttenergy[1]))]]), 0.] for ttenergy in datawf]
		randomset = [random.randint(0,len(datawf[0][1])-1) for i in range(len(datawf[0][1]))]
		wilson_flow_data_boot = [[ttenergy[0], numpy.mean([ttenergy[1][j] for j in randomset]), 0.] for ttenergy in datawf]

		wfdata_boot = [[(wilson_flow_data[t][0],),wilson_flow_data_boot[t][1],wilson_flow_data[t][2]] for t in range(td,tf)]
		
		minwf = opt.minimum_chisq(thWf,wfdata_boot,origwf,derivative_function=d_thWf)
		
		
		
		#[m, q, t_r, w_r]
		boot_data.append([minwf[0], minwf[1], (reference - minwf[1])/minwf[0], math.sqrt((reference)/minwf[0]), auto_t0(wilson_flow_data_boot), auto_w0(wilson_flow_data_boot)])
	return [[2*orig[i] - numpy.mean(zip(*boot_data)[i]),math.sqrt(float(len(datawf[0])-1)/len(datawf[0]))*numpy.std(zip(*boot_data)[i])] for i in range(6)]
	

wilson_fit = bootstrap_scale([[ttenergy[0][0], boot.block(zip(*ttenergy)[1][init:end],blocking)] for ttenergy in toplotwf],number_boot_resamplings)

title = "Wilson flow scale for run "+basename

import utils

resultstr = "Scale t_"+str(reference)+": "+utils.str_round(wilson_fit[2][0],wilson_fit[2][1])+" +/- "+utils.str_round(wilson_fit[2][1],wilson_fit[2][1])
print resultstr
title += "\n"+resultstr

resultstr = "Scale w_"+str(reference)+": "+utils.str_round(wilson_fit[3][0],wilson_fit[3][1])+" +/- "+utils.str_round(wilson_fit[3][1],wilson_fit[3][1])
print resultstr
title += "\n"+resultstr

resultstr = "Auto scale t_"+str(reference)+": "+utils.str_round(wilson_fit[4][0],wilson_fit[4][1])+" +/- "+utils.str_round(wilson_fit[4][1],wilson_fit[4][1])
print resultstr
title += "\n"+resultstr

resultstr = "Auto scale w_"+str(reference)+": "+utils.str_round(wilson_fit[5][0],wilson_fit[5][1])+" +/- "+utils.str_round(wilson_fit[5][1],wilson_fit[5][1])
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
	
coupling_str = ""
if measure_coupling:
	coupling = []
	for value in zip(wilson_flow_data, wilson_flow_data[1:]):
		if value[1][0] > tc:
			coupling = [((value[1][1]-value[0][1])/(value[1][0]-value[0][0]))*(tc-value[0][0])+value[0][1],((value[1][2]-value[0][2])/(value[1][0]-value[0][0]))*(tc-value[0][0])+value[0][2]]
			break
	print "Gauge coupling defined at ", tc, ": ", utils.str_round(coupling[0],coupling[1]), " +/- ", utils.str_round(coupling[1],coupling[1])
	coupling_str = ", "+utils.str_round(coupling[0],coupling[1])+", "+utils.str_round(coupling[1],coupling[1])
	
			

print "Mathematica output: {"+latticedata+(utils.str_round(wilson_fit[2][0],wilson_fit[2][1])+", "+utils.str_round(wilson_fit[2][1],wilson_fit[2][1])+", "+utils.str_round(wilson_fit[3][0],wilson_fit[3][1])+", "+utils.str_round(wilson_fit[3][1],wilson_fit[3][1])+", "+utils.str_round(wilson_fit[4][0],wilson_fit[4][1])+", "+utils.str_round(wilson_fit[4][1],wilson_fit[4][1])+", "+utils.str_round(wilson_fit[5][0],wilson_fit[5][1])+", "+utils.str_round(wilson_fit[5][1],wilson_fit[5][1])+coupling_str).replace("e","*10^")+"}"

try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.ticker as plticker

	plt.suptitle(title)

	ax = plt.subplot(2,1,1)
	ax.errorbar([elem[0] for elem in wilson_flow_data],[elem[1] for elem in wilson_flow_data],color='b',yerr=[elem[2] for elem in wilson_flow_data], marker='o',linestyle='')
	ax.plot([elem[0] for elem in wilson_flow_data],[thWf(elem[0],wilson_fit[0][0],wilson_fit[1][0]) for elem in wilson_flow_data],color='r')
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

	plt.subplots_adjust(hspace=0.3,wspace=0.3,left=0.2,top=0.8,bottom=0.1)
	plt.savefig(basename+"_wilson_flow")

except ImportError:
	print "Matplotlib not present, no plot output will be produced!"

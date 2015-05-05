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

def usage():
	print "Usage: pion_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: pion_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking --number_boot_resamplings=n --first_t_fit=t]"

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:b:s:k:", ["help", "folder=","basename=","skip=","blocking=","number_boot_resamplings=","first_t_fit="])
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
print "The error will be estimated with ", number_boot_resamplings, " bootstrapping resamplings (use --number_boot_resamplings to change it)"

import load_vector

data = load_vector.load(basename,"pion_exact.txt",folder)
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
import optimize as opt

propagator_data = [boot.bootstrap(boot.block(correlator[init:end],blocking)) for correlator in toplot]

Nt = len(propagator_data)

def thPm(t, T, m, A):
	return A*(math.exp(-m*t) + math.exp(m*(t - T)))

def d_thPm(t, T, m, A):
	return [A*(-t*math.exp(-m*t) + (t-T)*math.exp(m*(t - T))), (math.exp(-m*t) + math.exp(m*(t - T)))]

def chibt(m, A):
	return sum([((thPm(t, Nt, m, A) - propagator_data[t][0])**2)/(propagator_data[t][1]**2) for t in range(td-1,Nt-td)])

ppdata = [[(t,Nt),propagator_data[t][0],propagator_data[t][1]] for t in range(td-1,Nt-td)]

trial, minimum = opt.parallel_tempering(chibt, [1.,5.])
result = opt.minimum_chisq(thPm,ppdata,tuple(trial),derivative_function=d_thPm)
print "Pion mass: ",result[0]
print "Running bootstrapping analysis ..."
chisq = chibt(*result)/(Nt-td-3)

def bootstrap_correlator(data,trials=500):
	boot_data = []
	propagator_original = [boot.bootstrap(correlator) for correlator in data]
	ppdata_original = [[(t,Nt),propagator_original[t][0],propagator_original[t][1]] for t in range(td-1,Nt-td)]
	orig = opt.minimum_chisq(thPm,ppdata_original,result,derivative_function=d_thPm)
	for step in range(1,trials):
		propagator_boot = [numpy.mean([correlator[j] for j in [random.randint(0,len(correlator)-1) for i in range(len(correlator))]]) for correlator in data]
		ppdata_boot = [[(t,Nt),propagator_boot[t],propagator_original[t][1]] for t in range(td-1,Nt-td)]
		min = opt.minimum_chisq(thPm,ppdata_boot,orig,derivative_function=d_thPm)
		boot_data.append(list(min))
	return [[2*numpy.mean(orig[i]) - numpy.mean(zip(*boot_data)[i]),math.sqrt(float(len(data[0])-1)/len(data[0]))*numpy.std(zip(*boot_data)[i])] for i in range(2)]
	
		
pion_fit = bootstrap_correlator([boot.block(correlator[init:end],blocking) for correlator in toplot],number_boot_resamplings)

import utils

title = "Pion for run "+basename

roundlevel = int(-math.log10(pion_fit[0][1]) + 2)
resultstr = "Pion mass: "+utils.str_round(pion_fit[0][0],pion_fit[0][1])+" +/- "+utils.str_round(pion_fit[0][1],pion_fit[0][1])
print resultstr
title += "\n"+resultstr

resultstr = "Chi square: "+str(round(chisq,5))
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

print "Mathematica output: {"+latticedata+(utils.str_round(pion_fit[0][0],pion_fit[0][1])+", "+utils.str_round(pion_fit[0][1],pion_fit[0][1])).replace("e","*10^")+"}"

effective_mass = []
effective_mass_error = []
for i in range(len(propagator_data)-1):
	effective_mass.append(math.log(propagator_data[i][0]/propagator_data[i+1][0]))
	effective_mass_error.append((propagator_data[i + 1][0]/propagator_data[i][0])*math.sqrt(((propagator_data[i][1]/propagator_data[i + 1][0])**2) + (((propagator_data[i][0]/propagator_data[i + 1][0])**2)*((propagator_data[i + 1][1]/propagator_data[i + 1][0])**2))))

def effective_mass_function(x, mass):
		if x < Nt/2.:
			return mass
		else:
			return -mass
			
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.ticker as plticker
	
	plt.suptitle(title)

	ax = plt.subplot(2,1,1)
	ax.errorbar(range(len(effective_mass)),effective_mass,color='b',yerr=effective_mass_error, marker='o',linestyle='')
	ax.plot(numpy.arange(0,Nt/2,0.1),[effective_mass_function(x, pion_fit[0][0]) for x in numpy.arange(0,Nt/2,0.1)],color='r')
	ax.plot(numpy.arange(Nt/2,Nt,0.1),[effective_mass_function(x, pion_fit[0][0]) for x in numpy.arange(Nt/2,Nt,0.1)],color='r')
	ax.yaxis.set_label_text("Effective mass")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.2f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, Nt)

	log_correlator = []
	log_correlator_error = []
	for i in range(len(propagator_data)):
		log_correlator.append(math.log(propagator_data[i][0]))
		log_correlator_error.append((propagator_data[i][1]/propagator_data[i][0]))

	ax = plt.subplot(2,1,2)
	ax.errorbar(range(len(log_correlator)),log_correlator,color='b',yerr=log_correlator_error, marker='o',linestyle='')
	ax.plot(numpy.arange(0,Nt,0.1),[math.log(thPm(x, Nt, pion_fit[0][0], pion_fit[1][0])) for x in numpy.arange(0,Nt,0.1)],color='r')
	ax.yaxis.set_label_text("Log(correlator)")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.1f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, Nt)

	ax.xaxis.set_label_text("t")

	plt.subplots_adjust(wspace=0.3,left=0.2,top=0.85,bottom=0.1)
	plt.savefig(basename+"_pion")
except ImportError:
	print "Matplotlib not present, no plot output will be produced!"





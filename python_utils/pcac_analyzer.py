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
	print "Usage: pcac_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: pcac_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking --number_boot_resamplings=n --first_t_fit=t]"

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

datapp = load_vector.load(basename,"pion_exact.txt",folder)
finaldatapp = [j for i in datapp for j in i]

dataps = load_vector.load(basename,"ps_exact.txt",folder)
finaldataps = [j for i in dataps for j in i]

end=min(len(finaldataps),len(finaldatapp))

finaldatapp = finaldatapp[:end]
finaldataps = finaldataps[:end]

if end == 0:
	print "Data are not correctly loaded, do they exist?"
	sys.exit(2)
elif end < init:
	print "Not enough data to skip the first ",init," configs!"
	print "Use the option -s to a lower value than ", end
	sys.exit(2)
else:
	print "Analysing ", end-init, " configurations"

toplotps = zip(*finaldataps)
toplotpp = zip(*finaldatapp)

import bootstrapping as boot
import optimize as opt

propagator_pp_data = [boot.bootstrap(boot.block(correlator[init:end],blocking)) for correlator in toplotpp]
propagator_ps_data = [boot.bootstrap(boot.block(correlator[init:end],blocking)) for correlator in toplotps]

Nt = len(propagator_pp_data)

def thPm(t, T, m, A):
	return A*(math.exp(-m*t) + math.exp(m*(t - T)))

def d_thPm(t, T, m, A):
	return [A*(-t*math.exp(-m*t) + (t-T)*math.exp(m*(t - T))), (math.exp(-m*t) + math.exp(m*(t - T)))]

def chibtPm(m, A):
	return sum([((thPm(t, Nt, m, A) - propagator_pp_data[t][0])**2)/(propagator_pp_data[t][1]**2) for t in range(td-1,Nt-td)])

def thPs(t, T, m, A):
	return A*(math.exp(-m*t) - math.exp(m*(t - T)))

def d_thPs(t, T, m, A):
	return [A*(-t*math.exp(-m*t) - (t-T)*math.exp(m*(t - T))), (math.exp(-m*t) - math.exp(m*(t - T)))]

def chibtPs(m, A):
	return sum([((thPs(t, Nt, m, A) - propagator_ps_data[t][0])**2)/(propagator_ps_data[t][1]**2) for t in range(td-1,Nt-td)])

def thPc(t, T, m):
	return m

def d_thPc(t, T, m):
	return [ 1 ]


ppdata = [[(t,Nt),propagator_pp_data[t][0],propagator_pp_data[t][1]] for t in range(td-1,Nt-td)]
psdata = [[(t,Nt),propagator_ps_data[t][0],propagator_ps_data[t][1]] for t in range(td-1,Nt-td)]

trial, minimum = opt.parallel_tempering(chibtPm, [1.,5.])
resultpp = opt.minimum_chisq(thPm,ppdata,tuple(trial),derivative_function=d_thPm)

trial, minimum = opt.parallel_tempering(chibtPs, [1.,0.5])
resultps = opt.minimum_chisq(thPs,psdata,tuple(trial),derivative_function=d_thPs)

print "Pion mass: ",resultpp[0]
print "Running bootstrapping analysis ..."
chisq = chibtPm(*resultpp)/(Nt-td-3)
	
def bootstrap_pcac(datapp,dataps,trials=500):
	boot_data_pcac = []
	boot_data_pion = []
	boot_data_pcac_sym = []
	boot_data_gps = []
	boot_data_fps = []
	
	boot_data = []
	
	propagator_pp_original = [boot.bootstrap(correlator) for correlator in datapp]
	ppdata_original = [[(t,Nt),propagator_pp_original[t][0],propagator_pp_original[t][1]] for t in range(td-1,Nt-td)]
	propagator_ps_original = [boot.bootstrap(correlator) for correlator in dataps]
	psdata_original = [[(t,Nt),propagator_ps_original[t][0],propagator_ps_original[t][1]] for t in range(td-1,Nt-td)]
	pcacdata_original = [[(t,Nt),(propagator_ps_original[t+1][0] - propagator_ps_original[t-1][0])/(8.*propagator_pp_original[t][0]),math.sqrt((propagator_ps_original[t-1][1]**2 + propagator_ps_original[t+1][1]**2)/((8.*propagator_pp_original[t][0])**2) + ((propagator_pp_original[t][1]/(8.*propagator_pp_original[t][0]))*((propagator_ps_original[t+1][0] - propagator_ps_original[t-1][0])/(8.*propagator_pp_original[t][0])))**2)] for t in range(td-1,Nt-td)]
	
	origpp = opt.minimum_chisq(thPm,ppdata_original,resultpp,derivative_function=d_thPm)
	origps = opt.minimum_chisq(thPs,psdata_original,resultps,derivative_function=d_thPs)
	
	orig_pcac = origpp[0]*origps[1]/(4.*origpp[1])
	orig_pion = origpp[0]
	
	origpcac = opt.minimum_chisq(thPc,pcacdata_original,(orig_pcac,),derivative_function=d_thPc)
	orig_pcac_sym = origpcac[0]
	
	orig_gps = math.sqrt(origpp[0]*origpp[1])
	orig_fps = -origpcac[0]*math.sqrt(origpp[0]*origpp[1])/(origpp[0]**2)
	
	orig = [orig_pcac, -orig_pcac_sym, orig_pion, orig_gps, orig_fps]
	for step in range(1,trials):
		propagator_ps_boot = [boot.bootstrap([correlator[j] for j in [random.randint(0,len(correlator)-1) for i in range(len(correlator))]]) for correlator in dataps]
		propagator_pp_boot = [boot.bootstrap([correlator[j] for j in [random.randint(0,len(correlator)-1) for i in range(len(correlator))]]) for correlator in datapp]
		ppdata_boot = [[(t,Nt),propagator_pp_boot[t][0],propagator_pp_original[t][1]] for t in range(td-1,Nt-td)]
		psdata_boot = [[(t,Nt),propagator_ps_boot[t][0],propagator_ps_original[t][1]] for t in range(td-1,Nt-td)]
		minpp = opt.minimum_chisq(thPm,ppdata_boot,origpp,derivative_function=d_thPm)
		minps = opt.minimum_chisq(thPs,psdata_boot,origps,derivative_function=d_thPs)
		boot_data_pcac.append(minpp[0]*minps[1]/(4.*minpp[1]))
		
		boot_data_pion.append(minpp[0])
		
		
		pcacdata_boot = [[(t,Nt),(propagator_ps_boot[t+1][0] - propagator_ps_boot[t-1][0])/(8.*propagator_pp_boot[t][0]),pcacdata_original[t-td+1][1]] for t in range(td-1,Nt-td)]
		minpcac = opt.minimum_chisq(thPc,pcacdata_boot,(orig_pcac,),derivative_function=d_thPc)
		boot_data_pcac_sym.append(minpcac[0])
		
		boot_data_gps.append(math.sqrt(minpp[0]*minpp[1]))
		
		boot_data_fps.append(minpcac[0]*math.sqrt(minpp[0]*minpp[1])/(minpp[0]**2))
		
		#[pcac_mass_fit, pcac_sym_der, pion_mass, G_PS, F_PS]
		boot_data.append([minpp[0]*minps[1]/(4.*minpp[1]), -minpcac[0], minpp[0], math.sqrt(minpp[0]*minpp[1]), -(minpp[0]*minps[1]/(4.*minpp[1]))*math.sqrt(minpp[0]*minpp[1])/(minpp[0]**2)])
	return [[2*orig[i] - numpy.mean(zip(*boot_data)[i]),math.sqrt(float(len(datapp[0])-1)/len(datapp[0]))*numpy.std(zip(*boot_data)[i])] for i in range(5)]
	

pcac_fit = bootstrap_pcac([boot.block(correlator[init:end],blocking) for correlator in toplotpp],[boot.block(correlator[init:end],blocking) for correlator in toplotps],number_boot_resamplings)

#print pcac_fit

#pion_fit = bootstrap_correlator([boot.block(correlator[init:end],blocking) for correlator in toplotpp],number_boot_resamplings)

title = "Pion and PCAC mass for run "+basename

roundlevel = int(-math.log10(pcac_fit[2][1]) + 2)
resultstr = "Pion mass: "+str(round(pcac_fit[2][0],roundlevel))+" +/- "+str(round(pcac_fit[2][1],roundlevel))
print resultstr
title += "\n"+resultstr

resultstr = "Chi square (of pion mass fit): "+str(round(chisq,roundlevel))
print resultstr
title += "\n"+resultstr

roundlevel = int(-math.log10(pcac_fit[0][1]) + 2)
resultstr = "PCAC mass (fit definition): "+str(round(pcac_fit[0][0],roundlevel))+" +/- "+str(round(pcac_fit[0][1],roundlevel))
print resultstr
title += "\n"+resultstr

roundlevel = int(-math.log10(pcac_fit[1][1]) + 2)
resultstr = "PCAC mass (symmetric derivative definition): "+str(round(pcac_fit[1][0],roundlevel))+" +/- "+str(round(pcac_fit[1][1],roundlevel))
print resultstr
title += "\n"+resultstr

roundlevel = int(-math.log10(pcac_fit[3][1]) + 2)
resultstr = "G_PS: "+str(round(pcac_fit[3][0],roundlevel))+" +/- "+str(round(pcac_fit[3][1],roundlevel))
print resultstr
title += "\n"+resultstr

roundlevel = int(-math.log10(pcac_fit[4][1]) + 2)
resultstr = "Pion decay constant: "+str(round(pcac_fit[4][0],roundlevel))+" +/- "+str(round(pcac_fit[4][1],roundlevel))
print resultstr
title += "\n"+resultstr

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

	effective_mass = []
	effective_mass_error = []
	for i in range(len(propagator_pp_data)-1):
		effective_mass.append(math.log(propagator_pp_data[i][0]/propagator_pp_data[i+1][0]))
		effective_mass_error.append((propagator_pp_data[i + 1][0]/propagator_pp_data[i][0])*math.sqrt(((propagator_pp_data[i][1]/propagator_pp_data[i + 1][0])**2) + (((propagator_pp_data[i][0]/propagator_pp_data[i + 1][0])**2)*((propagator_pp_data[i + 1][1]/propagator_pp_data[i + 1][0])**2))))

	ax = plt.subplot(2,2,1)
	ax.errorbar(range(len(effective_mass)),effective_mass,color='b',yerr=effective_mass_error, marker='o',linestyle='')
	ax.plot(numpy.arange(0,Nt/2,0.1),[effective_mass_function(x, pcac_fit[2][0]) for x in numpy.arange(0,Nt/2,0.1)],color='r')
	ax.plot(numpy.arange(Nt/2,Nt,0.1),[effective_mass_function(x, pcac_fit[2][0]) for x in numpy.arange(Nt/2,Nt,0.1)],color='r')
	ax.yaxis.set_label_text("Effective pion\nmass")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.2f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, Nt)
	
	effective_mass = []
	effective_mass_error = []
	for t in range(2,len(propagator_pp_data)-1):
		effective_mass.append(-(propagator_ps_data[t+1][0] - propagator_ps_data[t-1][0])/(8.*propagator_pp_data[t][0]))
		effective_mass_error.append(math.sqrt((propagator_ps_data[t-1][1]**2 + propagator_ps_data[t+1][1]**2)/((8.*propagator_pp_data[t][0])**2) + ((propagator_pp_data[t][1]/(8.*propagator_pp_data[t][0]))*(-(propagator_ps_data[t+1][0] - propagator_ps_data[t-1][0])/(8.*propagator_pp_data[t][0])))**2))
	
	ax = plt.subplot(2,2,2)
	ax.errorbar(range(2,len(propagator_pp_data)-1),effective_mass,color='b',yerr=effective_mass_error, marker='o',linestyle='')
	ax.plot(numpy.arange(0,Nt,0.1),[pcac_fit[1][0] for x in numpy.arange(0,Nt,0.1)],color='r')
	ax.yaxis.set_label_text("PCAC mass")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.2f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, Nt)

	log_correlator = []
	log_correlator_error = []
	for i in range(len(propagator_pp_data)):
		log_correlator.append(math.log(propagator_pp_data[i][0]))
		log_correlator_error.append((propagator_pp_data[i][1]/propagator_pp_data[i][0]))

	ax = plt.subplot(2,2,3)
	ax.errorbar(range(len(log_correlator)),log_correlator,color='b',yerr=log_correlator_error, marker='o',linestyle='')
	ax.plot(numpy.arange(0,Nt,0.1),[math.log(thPm(x, Nt, pcac_fit[2][0], resultpp[1])) for x in numpy.arange(0,Nt,0.1)],color='r')
	ax.yaxis.set_label_text("Log(pion correlator)")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.1f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, Nt)
	
	ax.xaxis.set_label_text("t")
	
	correlator = []
	correlator_error = []
	for i in range(len(propagator_ps_data)):
		correlator.append(propagator_ps_data[i][0])
		correlator_error.append(propagator_ps_data[i][1])

	ax = plt.subplot(2,2,4)
	ax.errorbar(range(len(correlator)),correlator,color='b',yerr=correlator_error, marker='o',linestyle='')
	ax.plot(numpy.arange(0,Nt,0.1),[thPs(x, Nt, resultps[0], resultps[1]) for x in numpy.arange(0,Nt,0.1)],color='r')
	ax.yaxis.set_label_text("PS correlator")
	ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
	ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.1f'))
	ax.xaxis.set_major_locator(plticker.LinearLocator(5))
	ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.0f'))
	ax.set_xlim(0, Nt)

	ax.xaxis.set_label_text("t")

	plt.subplots_adjust(wspace=0.3,left=0.2,top=0.7,bottom=0.1)
	plt.savefig(basename+"_pcac")

except ImportError:
        print "Matplotlib not present, no plot output will be produced!"

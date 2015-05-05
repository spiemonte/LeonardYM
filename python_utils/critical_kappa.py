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

number_mcs = 100
filename=""

def usage():
	print "Usage: critical_kappa.py "
	print "Usage: critical_kappa.py --datafile=file [--number_mcs=n]"
	print "In file \"filename\", list the pion masses and their errors in the form [[kappa,pion_mass,error], ...]"

try:
	opts, args = getopt.getopt(sys.argv[1:], "", ["help", "number_mcs=", "datafile="])
except getopt.GetoptError, err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)
for option, argument in opts:
	if option in ("-h", "--help"):
		usage()
	elif option == "--number_mcs":
		number_mcs = int(argument)
	elif option == "--datafile":
		filename = argument

if filename == "":
	print "Filename not set!"
	usage()
	sys.exit(2)
	
print "The error will be estimated with ", number_mcs, " pseudo Monte Carlo samplings (use --number_mcs to change it)"

inputfile = open(filename,'r')
pion_masses=eval(inputfile.read())

pion_masses_sq = [[pion_mass[0], pion_mass[1]**2, 2.*pion_mass[1]*pion_mass[2]] for pion_mass in pion_masses]



if len(pion_masses) == 0:
	print "Data are not correctly loaded, try again!"
	sys.exit(2)

import optimize as opt

def thPm(x, m, c):
	return m*x+c

def d_thPm(x, m, c):
	return [x,1]

def chibt(m, c):
	return sum([((thPm(pion_mass_sq[0], m, c) - pion_mass_sq[1])**2)/(pion_mass_sq[2]**2) for pion_mass_sq in pion_masses_sq])

ppdata = [[(pion_mass_sq[0],),pion_mass_sq[1],pion_mass_sq[2]] for pion_mass_sq in pion_masses_sq]

trial, minimum = opt.parallel_tempering(chibt, [-1.,3.])
result = opt.minimum_chisq(thPm,ppdata,tuple(trial),derivative_function=d_thPm)
print "Critical kappa: ",result[1]
print "Running pseudo Monte Carlo analysis ..."
chisq = chibt(*result)/(len(pion_masses)-2)

def linear_pseudo_monte_carlo(data,trials=500):
	resampled_data = []
	for step in range(1,trials):
		ppdata_boot = [[(pion_mass_sq[0],),random.gauss(pion_mass_sq[1], pion_mass_sq[2]),pion_mass_sq[2]] for pion_mass_sq in data]
		min = opt.minimum_chisq(thPm,ppdata_boot,result,derivative_function=d_thPm)
		resampled_data.append(list(min))
	return [[numpy.mean(zip(*resampled_data)[i]),numpy.std(zip(*resampled_data)[i])] for i in range(2)]
	
		
critical_kappa_fit = linear_pseudo_monte_carlo(pion_masses_sq,number_mcs)

title = "Critical kappa fit"

roundlevel = int(-math.log10(critical_kappa_fit[0][1]) + 2)
resultstr = "Function: ("+str(round(critical_kappa_fit[0][0],roundlevel))+" +/- "+str(round(critical_kappa_fit[0][1],roundlevel))+")*x +"
roundlevel = int(-math.log10(critical_kappa_fit[1][1]) + 2)
resultstr += "("+str(round(critical_kappa_fit[1][0],roundlevel))+" +/- "+str(round(critical_kappa_fit[1][1],roundlevel))+")"
print resultstr
title += "\n"+resultstr

error=math.sqrt((critical_kappa_fit[1][1]/critical_kappa_fit[0][0])**2 + ((critical_kappa_fit[1][0]/critical_kappa_fit[0][0])**2)*((critical_kappa_fit[0][1]/critical_kappa_fit[0][0])**2))
roundlevel = int(-math.log10(error) + 2)
resultstr = "Critical kappa: "+str(round(-critical_kappa_fit[1][0]/critical_kappa_fit[0][0],roundlevel))+" +/- "+str(round(error,roundlevel))
print resultstr
title += "\n"+resultstr

resultstr = "Chi square: "+str(round(chisq,roundlevel))
print resultstr
title += "\n"+resultstr

plt.suptitle(title)

ax = plt.subplot(1,1,1)
ax.errorbar([pion_mass_sq[0] for pion_mass_sq in pion_masses_sq],[pion_mass_sq[1] for pion_mass_sq in pion_masses_sq],color='b',yerr=[pion_mass_sq[2] for pion_mass_sq in pion_masses_sq], marker='o',linestyle='')
xrange = numpy.arange(pion_masses_sq[0][0],pion_masses_sq[-1][0],0.00005)
ax.plot(xrange,[thPm(x, critical_kappa_fit[0][0],critical_kappa_fit[1][0]) for x in xrange],color='r')
ax.yaxis.set_label_text("(Pion mass)^2")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.2f'))
ax.xaxis.set_major_locator(plticker.LinearLocator(5))
ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.set_xlim(pion_masses_sq[0][0]-0.001,pion_masses_sq[-1][0]+0.001)

ax.xaxis.set_label_text("kappa")

plt.subplots_adjust(wspace=0.3,left=0.2,top=0.83,bottom=0.1)
plt.savefig("critical_kappa")






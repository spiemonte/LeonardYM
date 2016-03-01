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
gauge_factor = 0.
fermion_factor = 0.
subtraction_fermion = 0.
subtraction_gauge = 0.
error_subtraction_gauge = 0.
error_subtraction_fermion = 0.


def usage():
	print "Usage: trace_energy_momentum_analyzer.py -f folder -b basename [-s skip -k blocking]"
	print "Usage: trace_energy_momentum_analyzer.py --folder=folder --basename=basename [--skip=skip --blocking=blocking]"

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:b:s:k:::::::", ["help", "folder=","basename=","skip=","blocking=","fermion_factor=","gauge_factor=","subtraction_fermion=","subtraction_gauge=","error_subtraction_fermion=","error_subtraction_gauge="])
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
	elif option == "--fermion_factor":
		fermion_factor = float(argument)
	elif option == "--gauge_factor":
		gauge_factor = float(argument)
	elif option == "--subtraction_gauge":
		subtraction_gauge = float(argument)
	elif option == "--subtraction_fermion":
		subtraction_fermion = float(argument)
	elif option == "--error_subtraction_gauge":
		error_subtraction_gauge = float(argument)
	elif option == "--error_subtraction_fermion":
		error_subtraction_fermion = float(argument)
	
if basename == "":
	print "Basename not set!"
	usage()
	sys.exit(2)
if folder == "":
	print "Folder not set!"
	usage()
	sys.exit(2)
if gauge_factor == 0.:
	print "Use --gauge_factor= to set the inverse derivative of the logarithm beta function!"
	sys.exit(2)
if fermion_factor == 0.:
	print "Use --fermion_factor= to set the fermion mass with respect to \\beta!"
	sys.exit(2)
if subtraction_fermion == 0.:
	print "Use --subtraction_fermion= to set the chiral condensate at T=0!"
	sys.exit(2)
if subtraction_gauge == 0.:
	print "Use --subtraction_gauge= to set the gauge energy density at T=0!"
	sys.exit(2)
if error_subtraction_fermion == 0.:
	print "Warning, use --error_subtraction_fermion= to set the error on the chiral condensate at T=0!"
if error_subtraction_gauge == 0.:
	print "Warning, use --error_subtraction_gauge= to set the error on the gauge energy density at T=0!"
	
print "Skipping "+str(init)+" initial configurations as thermalization"
print "Using "+str(blocking)+" blocking levels"

import load_vector

datag = load_vector.load(basename,"gauge_energy.txt",folder)
datac = load_vector.load(basename,"condensate.txt",folder)

finaldatac = [j for i in datac for j in i]
finaldatag = [j for i in datag for j in i]

if len(finaldatac) < init:
	print "Not enough configurations the first ", init
	sys.exit(2)

finaldatag_cleaned = []
finaldatac_cleaned = []
removed = 0

for i in range(init,len(finaldatac)):
	if finaldatac[i][1] > 0.1:
		removed = removed + 1
	else:
		finaldatag_cleaned.append(finaldatag[i])
		finaldatac_cleaned.append(finaldatac[i])
		
finaldatag = finaldatag_cleaned
finaldatac = finaldatac_cleaned

print "Analyzing", len(finaldatag), "configurations"

print "Removing ", removed, " exceptional configurations."

import bootstrapping as boot
from utils import *

def bootstrap_trace_energy_momentum_tensor(datag,datac,volume,T,kappa,trials=500):
	boot_trace = []
	for i in range(trials):
		randomset = [random.randint(0,len(datag)-1) for i in range(1,len(datag))]
	
		subtraction_gauge_r = random.gauss(subtraction_gauge, error_subtraction_gauge)
		subtraction_fermion_r = random.gauss(subtraction_fermion, error_subtraction_fermion)
	
		gauge_energy_density_boot = numpy.mean([datag[j] for j in randomset])/volume
		chiral_condensate_boot = numpy.mean([datac[j] for j in randomset])

		boot_trace.append([-gauge_factor*(gauge_energy_density_boot - subtraction_gauge_r)*T*T*T*T + 2.*kappa*fermion_factor*gauge_factor*(chiral_condensate_boot - subtraction_fermion_r)*T*T*T*T/2., -gauge_factor*(gauge_energy_density_boot - subtraction_gauge_r)*T*T*T*T, - 2.*kappa*fermion_factor*gauge_factor*(chiral_condensate_boot - subtraction_fermion_r)*T*T*T*T/2.]);
	return [[numpy.mean(zip(*boot_trace)[i]), numpy.std(zip(*boot_trace)[i])] for i in range(3)]

volume = 0
T = 0
kappa = 0.
pattern_volume = re.search('([0-9]+)c([0-9]+)',basename)
if pattern_volume != None:
	volume = float(int(pattern_volume.group(1))*int(pattern_volume.group(1))*int(pattern_volume.group(1))*int(pattern_volume.group(2)))
	T = float(int(pattern_volume.group(2)))
pattern_kappa = re.search('([0-9]+)k',basename)
if pattern_kappa != None:
	kappa = float("0."+pattern_kappa.group(1))

energy_momentum_tensor = bootstrap_trace_energy_momentum_tensor(boot.block(zip(*finaldatag)[0],blocking),boot.block(zip(*finaldatac)[0],blocking),volume,T,kappa)

resultstr = "Trace of the energy momentum tensor: "+str_round(energy_momentum_tensor[0][0],energy_momentum_tensor[0][1])+" +/- "+str_round(energy_momentum_tensor[0][1],energy_momentum_tensor[0][1])+"\n"
resultstr += "Trace of the energy momentum tensor (gauge contribution): "+str_round(energy_momentum_tensor[1][0],energy_momentum_tensor[1][1])+" +/- "+str_round(energy_momentum_tensor[1][1],energy_momentum_tensor[1][1])+"\n"
resultstr += "Trace of the energy momentum tensor (fermion contribution): "+str_round(energy_momentum_tensor[2][0],energy_momentum_tensor[2][1])+" +/- "+str_round(energy_momentum_tensor[2][1],energy_momentum_tensor[2][1])

print resultstr

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

print "Mathematica output: {"+latticedata+(str_round(energy_momentum_tensor[0][0],energy_momentum_tensor[0][1])+", "+str_round(energy_momentum_tensor[0][1],energy_momentum_tensor[0][1])+", "+str_round(energy_momentum_tensor[1][0],energy_momentum_tensor[1][1])+", "+str_round(energy_momentum_tensor[1][1],energy_momentum_tensor[1][1])+", "+str_round(energy_momentum_tensor[2][0],energy_momentum_tensor[2][1])+", "+str_round(energy_momentum_tensor[2][1],energy_momentum_tensor[2][1])).replace("e","*10^")+"}"








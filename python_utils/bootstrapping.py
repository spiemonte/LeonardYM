#!/usr/bin/python
import numpy
import random
import math

def block(data, recursion):
	for r in range(recursion):
		tmp = []
		for i in range(len(data)%2,len(data),2):
			tmp.append((data[i] + data[i+1])/2.)
		data = tmp
	return data

def bootstrap_susc(data, data_sq, number_resamplings=500):
	boot = []
	for i in range(1,number_resamplings):
		tmp = [[data[j], data_sq[j]] for j in [random.randint(0,len(data)-1) for i in range(1,len(data))]]
		tmp = zip(*tmp)
		boot.append(numpy.mean(tmp[1]) - numpy.mean(tmp[0])*numpy.mean(tmp[0]))
	return [2*(numpy.mean(data_sq) - numpy.mean(data)*numpy.mean(data)) - numpy.mean(boot),math.sqrt(float(len(data)-1)/len(data))*numpy.std(boot)]

def bootstrap(data, number_resamplings=500):
	boot = []
	for i in range(1,number_resamplings):
		boot.append(numpy.mean([data[j] for j in [random.randint(0,len(data)-1) for i in range(1,len(data))]]))
	return [2*numpy.mean(data) - numpy.mean(boot),math.sqrt(float(len(data)-1)/len(data))*numpy.std(boot)]

version = '1.1'

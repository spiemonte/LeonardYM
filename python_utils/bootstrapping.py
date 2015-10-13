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

#def bootstrap_binder(data, data_sq, data_th, data_ft, number_resamplings=500):
#	boot = []
#	for i in range(1,number_resamplings):
#		tmp = [[data[j], data_sq[j], data_th[j], data_ft[j]] for j in [random.randint(0,len(data)-1) for i in range(1,len(data))]]
#		tmp = zip(*tmp)
#		mu = numpy.mean(tmp[0])
#		sq = numpy.mean(tmp[1])
#		th = numpy.mean(tmp[2])
#		ft = numpy.mean(tmp[3])
#		boot.append(1. - (ft + 6.*(mu**2)*sq - 4.*mu*th-3.*(mu**4))/(3.*((sq - mu**2)**2)) )
#	mu = numpy.mean(data)
#	sq = numpy.mean(data_sq)
#	th = numpy.mean(data_th)
#	ft = numpy.mean(data_ft)
#	print "Original", (1. - (ft + 6.*(mu**2)*sq - 4.*mu*th-3.*(mu**4))/(3.*((sq - mu**2)**2)))
#	print "Derived", numpy.mean(boot)
#	return [2.*(1. - (ft + 6.*(mu**2)*sq - 4.*mu*th-3.*(mu**4))/(3.*((sq - mu**2)**2))) - numpy.mean(boot),math.sqrt(float(len(data)-1)/len(data))*numpy.std(boot)]

def bootstrap_binder(data_sq, data_ft, number_resamplings=500):
	boot = []
	for i in range(1,number_resamplings):
		tmp = [[data_sq[j], data_ft[j]] for j in [random.randint(0,len(data_sq)-1) for i in range(1,len(data_sq))]]
		tmp = zip(*tmp)
		sq = numpy.mean(tmp[0])
		ft = numpy.mean(tmp[1])
		boot.append(1. - (ft)/(3.*sq*sq) )
	sq = numpy.mean(data_sq)
	ft = numpy.mean(data_ft)
	return [2.*(1. - (ft)/(3.*sq*sq)) - numpy.mean(boot),math.sqrt(float(len(data_sq)-1)/len(data_sq))*numpy.std(boot)]

def bootstrap(data, number_resamplings=500):
	boot = []
	for i in range(1,number_resamplings):
		boot.append(numpy.mean([data[j] for j in [random.randint(0,len(data)-1) for i in range(1,len(data))]]))
	return [2*numpy.mean(data) - numpy.mean(boot),math.sqrt(float(len(data)-1)/len(data))*numpy.std(boot)]

version = '1.1'

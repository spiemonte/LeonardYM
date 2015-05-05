#!/usr/bin/python
import random
import math
import numpy

def derivative(f,x,h=0.01):
	a = [[0 for i in range(10)] for j in range(10)]
	a[0][0] = (f(x+h)-f(x-h))/(2*h)
	error = 1000000.
	answer = 0.
	for i in range(10):
		h = h/1.4
		a[0][i] = (f(x + h) - f(x - h)) / (2 * h);
		fac = 1.4 * 1.4;
		for j in range(1,i+1):
			a[j][i]=(a[j-1][i] * fac - a[j-1][i-1]) / (fac - 1)
			fac *= 1.4 * 1.4
			err = max(abs(a[j][i] - a[j-1][i]), abs(a[j][i] - a[j-1][i-1]))
			if (err <= error):
				error = err;
				answer = a[j][i];
		if (abs(a[i][i] - a[i-1][i-1]) >= 2 * error):
			break
	return answer

def gradient(f,point):
	gradient = []
	for i in range(len(point)):
		ftmp = lambda x: f(point[0:i]+[x]+point[i+1:len(point)])
		gradient.append(derivative(ftmp,point[i]))
	return gradient

def minimum_chisq(function, points, initial_guess,derivative_function=None,tolerance=0.0000000001,maxiterations=1000):
	#We use the Levenberg-Marquardt algorithm
	result = initial_guess
	for step in range(maxiterations):
		J = []
		for point in points:
			if derivative_function == None:
				row = []
				for i in range(len(result)):
					ftmp = lambda x: function(*(point[0]+result[0:i]+(x,)+result[i+1:len(result)]))
					row.append(derivative(ftmp,result[i])/point[2])
				J.append(row)
			else:
				d = derivative_function(*(point[0]+result));
				J.append([elem/point[2] for elem in d])
		R = [(function(*(point[0]+result)) - point[1])/point[2] for point in points]
		JM = numpy.array(J).transpose()
		RM = numpy.array(R)
		V = numpy.dot(JM,RM)
		A = numpy.dot(JM,JM.transpose())
		solution = numpy.linalg.solve(A, V)
		next = [result[i] - solution[i] for i in range(len(solution))]
		result = tuple(next)
		norm = math.sqrt(sum([elem*elem for elem in solution]))
		if norm < tolerance:
			break
		elif step == maxiterations-1:
			print "Warning, failure in finding convergence for minimum_chisq!"
	return result

def metropolis(delta):
	if delta < -100.:
		return False
	elif delta >= 0.:
		return True
	elif math.exp(delta) > random.uniform(0,1):
		return True
	else:
		return False

def parallel_tempering(function, initial_guess, beta = [5.,1.,0.2], epsilon = [0.003,0.01,0.1], steps = 2000):
	points = [[x + random.uniform(-1,1)*epsilon[i] for x in initial_guess] for i in range(len(beta))]
	fpoints = [function(*tuple(point)) for point in points]
	result = initial_guess
	minimum = 100000000000
	acceptance = 0
	swaps = 0
	for step in range(steps):
		for i in range(len(epsilon)):
			#metropolis search
			move = [random.uniform(-1,1)*epsilon[i] for j in range(len(points[i]))]
			trial = [points[i][j] + move[j] for j in range(len(points[i]))]
			ftrial = function(*tuple(trial))
			if metropolis(beta[i]*(fpoints[i]-ftrial)):
				points[i] = trial
				fpoints[i] = ftrial
				acceptance = acceptance + 1
				if fpoints[i] < minimum:
					minimum = fpoints[i]
					result = points[i]
		#parallel tempering
		for i in range(len(epsilon)-1):
			if metropolis((beta[i]-beta[i+1])*(fpoints[i]-fpoints[i+1])):
				#Swap
				points[i], points[i+1] = points[i+1], points[i]
				fpoints[i], fpoints[i+1] = fpoints[i+1], fpoints[i]
				swaps = swaps + 1
	print "MC acceptance: "+str(float(acceptance)/(steps*len(epsilon)))
	print "PT swap rate: "+str(float(swaps)/(steps*(len(epsilon)-1)))
	return (result,minimum)


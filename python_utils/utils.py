#!/usr/bin/python

import math
import numpy

def str_round(number,error):
	return str(round(number, int(-math.log10(error) + 2)))

class Polynomial:
	def __init__(self, coefficients = [0]):
		self.coefficients = coefficients
	def __add__(self, toAdd):
		result = Polynomial();
		if len(self.coefficients) >= len(toAdd.coefficients):
			result.coefficients = self.coefficients
			for i in range(len(toAdd.coefficients)):
				result.coefficients[i] += toAdd.coefficients[i]
		else:
			result.coefficients = toAdd.coefficients
			for i in range(len(self.coefficients)):
				result.coefficients[i] += self.coefficients[i]
		return result
		
	def __mul__(self, toMul):
		result = Polynomial();
		for k in range(len(self.coefficients)):
			ctmp = [0 for i in range(k)]
			for m in range(len(toMul.coefficients)):
				ctmp.append(toMul.coefficients[m]*self.coefficients[k])
			result = result + Polynomial(ctmp)
		return result
		
	def derive(self):
		ctmp = []
		for c in list(enumerate(self.coefficients))[1:]:
			ctmp.append(c[0]*c[1])
		return Polynomial(ctmp)
		
	def zero(self,a,b,niters=500):
		if self.evaluate(a) < 0 and self.evaluate(b) < 0:
			raise NameError('No zeros!')
		elif self.evaluate(a) > 0 and self.evaluate(b) > 0:
			raise NameError('No zeros!')
		
		x0 = a
		y0 = self.evaluate(a)
		
		x1 = b
		y1 = self.evaluate(b)
		
		for i in range(niters):
			x = x1 - (x1-x0)*y1/(y1-y0)
			y = self.evaluate(x);
			
			if abs(y) < 10e-13:
				break
			
			if y*y1 > 0:
				y1 = y
				x1 = x
				
			else:
				y0 = y
				x0 = x
				
		if abs(y1) < abs(y0):
			return x1
		else:
			return x0
		
	def evaluate(self, x):
		result = self.coefficients[-1]
		for c in reversed(self.coefficients[:-1]):
			result = result*x + c
		return result
		
	def __str__(self):
		res = "";
		for k in range(len(self.coefficients)):
			res = res+str(self.coefficients[k])+"x^"+str(k)+" +"
		return res[0:len(res)-2]
		
def LagrangePolynomial(points,values):
	result = Polynomial()
	for i in range(len(points)):
		v = 1.;
		for point in points[:i]+points[i+1:]:
			v = v*(points[i]-point)
		tmp = Polynomial([values[i]/v])
		for point in points[:i]+points[i+1:]:
			tmp = tmp*Polynomial([-point,1])
		result = result + tmp
	return result

def LeastSquarePolynomial(points,values,order):
	X = [[math.pow(p,j) for j in range(0,order+1)] for p in points]
	XM = numpy.array(X)
	YM = numpy.array(values)
	solution = numpy.linalg.solve(numpy.dot(XM.transpose(),XM), numpy.dot(XM.transpose(),YM))
	result = Polynomial([a for a in solution])
	return result
	
	

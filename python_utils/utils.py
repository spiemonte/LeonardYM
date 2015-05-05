#!/usr/bin/python

import math

def str_round(number,error):
	return str(round(number, int(-math.log10(error) + 2)))
	


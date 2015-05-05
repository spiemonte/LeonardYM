#!/usr/bin/python
import os
import sys

if len(sys.argv) != 4:
        print "This script must be called with two arguments!"
        print "Usage:\n   vector_load.py file data folder"
        sys.exit(1)

def load_vector(basename,endname,folder,max):
	result = []
	level = 0
	tmp = ""
	for index in range(1,max):
		if os.path.isfile(folder+basename+"_"+str(index)+"_"+endname):
			inputfile = open(folder+basename+"_"+str(index)+"_"+endname,'r')
			filedata = inputfile.read()
			for character in filedata:
				if character == "{":
					if level == 0:
						result.append([])
					if level == 1:
						result[-1].append([])
					if level == 2:
						result[-1][-1].append([])
					if level == 3:
						result[-1][-1][-1].append([])
					level = level + 1
				elif character == "}":
					if tmp != "":
						if level == 0:
							result.append(float(tmp))
						elif level == 1:
							result[-1].append(float(tmp))
						elif level == 2:
							result[-1][-1].append(float(tmp))
						elif level == 3:
							result[-1][-1][-1].append(float(tmp))
						tmp = ""
					level = level - 1
					if level < 0:
						level = 0
				elif character == ",":
					if tmp != "":
						if level == 0:
							result.append(float(tmp))
						elif level == 1:
							result[-1].append(float(tmp))
						elif level == 2:
							result[-1][-1].append(float(tmp))
						elif level == 3:
							result[-1][-1][-1].append(float(tmp))
						tmp = ""
				else:
					if unicode(character).isnumeric():
						tmp = tmp + character
					elif character == "-":
						tmp = tmp + character
					elif character == ".":
						tmp = tmp + character
					elif character == "e":
						tmp = tmp + character
			if tmp != "":
						if level == 0:
							result.append(float(tmp))
						elif level == 1:
							result[-1].append(float(tmp))
						elif level == 2:
							result[-1][-1].append(float(tmp))
						elif level == 3:
							result[-1][-1][-1].append(float(tmp))
						tmp = ""
	return result

	

data = 	load_vector(sys.argv[1],sys.argv[2],sys.argv[3],10)
print data
		
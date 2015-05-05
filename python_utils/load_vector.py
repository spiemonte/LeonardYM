#!/usr/bin/python
import os


def load(basename,endname,folder,log=False):
	result = []
	level = 0
	tmp = ""
	skipped = 0
	for index in range(1,50000):
		if os.path.isfile(folder+basename+"_"+str(index)+"_"+endname):
			inputfile = open(folder+basename+"_"+str(index)+"_"+endname,'r')
			filedata = inputfile.read()
			filedata.replace("},\n{","},")
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
			inputfile.close()
			if level != 1:
				print "Warning, possible corrupted data in",folder+basename+"_"+str(index)+"_"+endname
				level = 0
		else:
			if log:
				print "Skipping file "+folder+basename+"_"+str(index)+"_"+endname
			skipped = skipped + 1
		level = 0
		if skipped == 10000:
			break
	return result

version = '1.1'

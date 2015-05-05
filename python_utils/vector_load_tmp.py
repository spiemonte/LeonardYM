#!/usr/bin/python
import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker


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
		else:
			print "Skipping file "+folder+basename+"_"+str(index)+"_"+endname
		level = 0
		
	return result

	
def get_ticks(data,nticks,prec):
	return [round(2*(i+(max(data) - min(data)) / (2*nticks)),prec-1)/2 for i in np.arange(round(2*min(data),prec-1)/2, max(data), round(2*((max(data) - min(data)) / nticks),prec-1)/2)]	

data = 	load_vector(sys.argv[1],sys.argv[2],sys.argv[3],30)
finaldata = [j for i in data for j in i]

init=100
end=len(finaldata)

toplot = zip(*finaldata)
ax = plt.subplot(3,2,1)
ax.plot(range(init,end),toplot[0][init:end],color='b')
ax.yaxis.set_label_text("Total\nplaquette")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.xaxis.set_major_locator(plticker.LinearLocator(5))

ax = plt.subplot(3,2,3)
ax.plot(range(init,end),toplot[1][init:end],color='r')
ax.yaxis.set_label_text("Spatial\nplaquette")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.xaxis.set_major_locator(plticker.LinearLocator(5))

ax = plt.subplot(3,2,5)
ax.plot(range(init,end),toplot[2][init:end],color='g')
ax.yaxis.set_label_text("Temporal\nplaquette")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))
ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.xaxis.set_major_locator(plticker.LinearLocator(5))

ax.xaxis.set_label_text("T Monte Carlo")

ax = plt.subplot(3,2,2)
ax.hist(toplot[0][init:end],facecolor='b')
ax.xaxis.set_major_locator(plticker.LinearLocator(3))
ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.yaxis.set_label_text("Frequency")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

ax = plt.subplot(3,2,4)
ax.hist(toplot[1][init:end],facecolor='r')
ax.xaxis.set_major_locator(plticker.LinearLocator(3))
ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.yaxis.set_label_text("Frequency")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

ax = plt.subplot(3,2,6)
ax.hist(toplot[2][init:end],facecolor='g')
ax.xaxis.set_major_locator(plticker.LinearLocator(3))
ax.xaxis.set_major_formatter(plticker.FormatStrFormatter('%0.4f'))
ax.xaxis.set_label_text("Plaquette value")
ax.yaxis.set_label_text("Frequency")
ax.yaxis.set_major_locator(plticker.MaxNLocator(6,prune='both'))

plt.subplots_adjust(wspace=0.3,left=0.2)
plt.savefig('myfig')
		






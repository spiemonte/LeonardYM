nc = 2

CC       = gcc
CPP      = g++
OPTIMIZE = -O3 -march=native -Wall -Wextra -pipe -fopenmp
INC	= -I./extra_libs/Eigen -I./extra_libs/
CFLAGS   = $(OPTIMIZE) $(INC)
CPPFLAGS = $(OPTIMIZE) $(INC) -DMBE="\"$(PWD)/source/MatrixBaseExtension.h\"" -DNUMCOLORS=$(nc) -DMULTITHREADING -DADJOINT

include Makefile.obj.mk

all: leonardQCD
	$(CPP) $(CPPFLAGS) $(OBJECTS) -o leonardQCD.exe -lboost_program_options
	@echo "Compiled with NUMCOLORS: " $(nc)

leonardQCD: $(OBJECTS)
	@echo "done!"
	
include Makefile.list.mk

clean:
	-rm ./*.exe
	-rm ./*.o



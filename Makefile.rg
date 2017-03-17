nc = 2

ifeq ($(adjoint),false)
	adj = 
else
	adj = -DADJOINT
endif

CC       = gcc
CPP      = g++
OPTIMIZE = -O2 -Wno-unused-local-typedefs -Wall -Wextra -pipe -fopenmp -mssse3
INC	= -I./extra_libs/Eigen -I./extra_libs/ -I./source/ 
CFLAGS   = $(OPTIMIZE) $(INC)
CPPFLAGS = $(OPTIMIZE) $(INC) -DMBE="\"$(PWD)/source/utils/MatrixBaseExtension.h\"" -DNUMCOLORS=$(nc) $(adj) -DEIGEN -DMULTITHREADING

include Makefile.obj.mk

all: leonardQCD
	$(CPP) $(CPPFLAGS) $(OBJECTS) -o leonardQCD.exe /usr/lib/x86_64-linux-gnu/libboost_program_options.a 
	@echo "Compiled with NUMCOLORS: " $(nc)

leonardQCD: $(OBJECTS)
	@echo "done!"
	
include Makefile.list.mk

clean:
	-rm ./*.exe
	-rm ./*.o



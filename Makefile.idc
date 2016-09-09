nc = 2

ifeq ($(adjoint),false)
	adj = 
else
	adj = -DADJOINT
endif

CC       = gcc
CPP      = g++
OPTIMIZE = -O2 -Wno-unused-local-typedefs -Wall -Wextra -pipe -fopenmp -mssse3
INC	= -I./extra_libs/Eigen -I./extra_libs/ -I./source/ -I/home/pis52920/local/usr/include/
CFLAGS   = $(OPTIMIZE) $(INC)
CPPFLAGS = $(OPTIMIZE) $(INC) -DMBE="\"$(PWD)/source/utils/MatrixBaseExtension.h\"" -DNUMCOLORS=$(nc) $(adj) -DEIGEN -DMULTITHREADING

include Makefile.obj.mk

all: leonardQCD
	$(CPP) $(CPPFLAGS) $(OBJECTS) -o leonardQCD.exe /home/pis52920/local/usr/lib/libboost_program_options.so.1.49.0 -L/home/pis52920/local/usr/lib/ -Wl,-rpath,/home/pis52920/local/usr/lib/
	@echo "Compiled with NUMCOLORS: " $(nc)

leonardQCD: $(OBJECTS)
	@echo "done!"
	
include Makefile.list.mk

clean:
	-rm ./*.exe
	-rm ./*.o



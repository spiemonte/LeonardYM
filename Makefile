nc = 2

ifeq ($(adjoint),false)
	adj = 
else
	adj = -DADJOINT
endif

CC       = gcc
CPP      = g++-7
OPTIMIZE = -O2 -march=native -Wall -Wextra -pipe -fopenmp -Wno-int-in-bool-context -L/usr/local/boost-1.66.0/lib/
INC	= -I./extra_libs/Eigen -I./extra_libs/ -I./source/ -I/usr/local/boost-1.66.0/include/
CFLAGS   = $(OPTIMIZE) $(INC)
CPPFLAGS = $(OPTIMIZE) $(INC) -DMBE="\"$(PWD)/source/utils/MatrixBaseExtension.h\"" -DNUMCOLORS=$(nc) -DMULTITHREADING $(adj) -DEIGEN

include Makefile.obj.mk

all: leonardQCD
	$(CPP) $(CPPFLAGS) $(OBJECTS) -o leonardQCD.exe /usr/local/boost-1.66.0/lib/libboost_program_options.a
	@echo "Compiled with NUMCOLORS: " $(nc)

leonardQCD: $(OBJECTS)
	@echo "done!"
	
include Makefile.list.mk

clean:
	-rm ./*.exe
	-rm ./*.o



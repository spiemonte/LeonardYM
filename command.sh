#!/bin/sh

mpirun -n 8 -hosts mic0,mic1,mic2,mic3,mic4,mic5,mic6,mic7 ./leonardQCD.exe --configfile=config_nf32_multishift_three_levels_two_pseudofermions.cfg

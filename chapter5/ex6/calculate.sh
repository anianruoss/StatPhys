#!/bin/bash

# prerequisites:
# - current directory contains executables
#   * mdatom
# - gr.py (with execute permission)

DIR=calculations_dir
rm -rf $DIR
mkdir $DIR
cd $DIR

# Perform periodic calculation
./../mdatom ../params_periodic.inp ../coord.inp > periodic.out
# radial distribution
./../gr.py periodic.out > gr_periodic.dat

# Perform vacuum calculation, equilibration
./../mdatom ../params_vacuum_eq.inp ../coord.inp > vacuum_eq.out
cp coords.final coords_eq.inp

# vacuum calculation
./../mdatom ../params_vacuum.inp coords_eq.inp > vacuum.out
# radial distribution
./../gr.py vacuum.out > gr_vacuum.dat

cd ..

#!/bin/bash

# prerequisites:
# - current directory contains executables
#   * mdatom

DIR=calculations_dir
rm -rf $DIR
mkdir $DIR

cd $DIR
for TT in 0.003 0.01 0.03 0.1 0.3 1.0
do
	# Create parameter file
	PARAM_FILE=$TT.inp
	sed "s/_CouplingTime_/$TT/" < ../params.template > $PARAM_FILE

	# Perform calculations
	OUTPUT_FILE=${TT}.out
	./../mdatom $PARAM_FILE ../coord.inp > $OUTPUT_FILE
done
cd ..

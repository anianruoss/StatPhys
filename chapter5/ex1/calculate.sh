#!/bin/bash

# prerequisites:
# - current directory contains executables
#   * mdatom_double
#   * mdatom_float
#   * mdatom_long

DIR=calculations_dir
rm -rf $DIR
mkdir $DIR

cd $DIR
for STEPS in 100 200 500 1000 2000 5000 10000 20000 50000
do
	# Create parameter file
	PARAM_FILE=$STEPS.inp
	sed s/_NumberMDSteps_/$STEPS/ < ../params.template > $PARAM_FILE

	# Perform calculations for the three types
	for TYPE in double float long
	do
		OUTPUT_FILE=${TYPE}_${STEPS}.out
		./../mdatom_$TYPE $PARAM_FILE ../coord.inp > $OUTPUT_FILE
	done
done
cd ..

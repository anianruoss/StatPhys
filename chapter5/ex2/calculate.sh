#!/bin/bash

# prerequisites:
# - current directory contains executables
#   * mdatom

DIR=calculations_dir
rm -rf $DIR
mkdir $DIR

cd $DIR
for DT in 0.02 0.01 0.005 0.002 0.001
do
	# set number of steps so that total time is 20
	STEPS="$(echo "20/${DT}" | bc)"
	# Create parameter file
	PARAM_FILE=$DT.inp
	sed "s/_NumberMDSteps_/$STEPS/; s/_TimeStep_/$DT/" < ../params.template > $PARAM_FILE

	# Perform calculations
	OUTPUT_FILE=${DT}.out
	./../mdatom $PARAM_FILE ../coord.inp > $OUTPUT_FILE
done
cd ..

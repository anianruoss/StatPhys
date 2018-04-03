#!/bin/bash

# prerequisites:
# - current directory contains executables
#   * mdatom

DIR=calculations_dir
rm -rf $DIR
mkdir $DIR

cd $DIR
for RCUTF in 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0
do
	# Create parameter file
	PARAM_FILE=$RCUTF.inp
	sed "s/_Cutoff_/$RCUTF/" < ../params.template > $PARAM_FILE

	# Perform calculations
	OUTPUT_FILE=${RCUTF}.out
	./../mdatom $PARAM_FILE ../coord.inp > $OUTPUT_FILE
done
cd ..

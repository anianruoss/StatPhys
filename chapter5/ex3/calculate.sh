#!/bin/bash

# prerequisites:
# - current directory contains executables
#   * mdatom
#   * gr.py (must have execute rights)

DIR=calculations_dir
rm -rf $DIR
mkdir $DIR

cd $DIR
for TYPE in "constant_volume" "constant_density"
do
	mkdir $TYPE
	cd $TYPE
	for NCUBICROOT in 5 6 7 8 9 10
	do
		# Calculate total number of atoms
		NATOMSTOT="$(echo "${NCUBICROOT} * ${NCUBICROOT} * ${NCUBICROOT}" | bc)"
		# Calculate volume for constant density
		if [ $TYPE == "constant_density" ]; then
			BOXSIZE="$(echo "${NCUBICROOT}*2" | bc)"
		else
			BOXSIZE="10"
		fi
		# Create parameter file
		PARAM_FILE_EQ=${NATOMSTOT}_eq.inp
		PARAM_FILE_RUN=${NATOMSTOT}_run.inp
		sed "s/_NumberAtoms_/$NATOMSTOT/g; s/_BoxSize_/$BOXSIZE/g; s/_NAtomsEdge_/$NCUBICROOT/g" < ../../params_eq.template > $PARAM_FILE_EQ
		sed "s/_NumberAtoms_/$NATOMSTOT/g; s/_BoxSize_/$BOXSIZE/g; s/_NAtomsEdge_/$NCUBICROOT/g" < ../../params_run.template > $PARAM_FILE_RUN

		# Perform equilibration calculation
		OUTPUT_FILE_EQ=${NATOMSTOT}_eq.out
		./../../mdatom $PARAM_FILE_EQ > $OUTPUT_FILE_EQ
		cp coords.final coord_${NATOMSTOT}_eq.inp
		# Perform simulation
		OUTPUT_FILE_RUN=${NATOMSTOT}_run.out
		./../../mdatom $PARAM_FILE_RUN coord_${NATOMSTOT}_eq.inp > $OUTPUT_FILE_RUN

		# radial distribution
		./../../gr.py $OUTPUT_FILE_RUN > gr_${NATOMSTOT}.dat
	done
	cd ..
done
cd ..

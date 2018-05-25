#!/usr/bin/env bash

SOURCES=../../mdatom
BUILD=$SOURCES/build
UTIL=../../util

source ../../virtualenv/bin/activate

DIR=data
rm -rf $DIR
mkdir $DIR
cd $DIR

for TEMPERATURE in 0 100 200 300 400
do
    python ../$UTIL/pargen.py ../mdpar.ini > params_t$TEMPERATURE.inp InitialTemperature=$TEMPERATURE
done

for NUMATOMS in 2 4 8 16 32
do
    python ../$UTIL/pargen.py ../mdpar.ini > params_n$NUMATOMS.inp NumberAtoms=$NUMATOMS
done

for SIMTYPE in Base Harmonic Shake
do
    for TEMPERATURE in 0 100 200 300 400
    do
        ../$BUILD/mdatom$SIMTYPE params_t$TEMPERATURE.inp ../coords.inp > out_${SIMTYPE}_t$TEMPERATURE
        mv coords.traj coords_${SIMTYPE}_t$TEMPERATURE.traj
    done
    for NUMATOMS in 2 4 8 16 32
    do
        ../$BUILD/mdatom$SIMTYPE params_n$NUMATOMS.inp ../coords.inp > out_${SIMTYPE}_n$NUMATOMS
        mv coords.traj coords_${SIMTYPE}_n$NUMATOMS.traj
    done
done

cd ..

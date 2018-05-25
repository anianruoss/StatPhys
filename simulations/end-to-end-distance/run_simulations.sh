#!/usr/bin/env bash

SOURCES=../../mdatom
BUILD=$SOURCES/build
UTIL=../../util

source ../../virtualenv/bin/activate

DIR=calculations_dir
rm -rf $DIR
mkdir $DIR
cd $DIR

for TEMPERATURE in 0 100 200 300 400
do
    python ../$UTIL/pargen.py ../mdpar.ini > params_t$TEMPERATURE.inp InitialTemperature=$TEMPERATURE
done

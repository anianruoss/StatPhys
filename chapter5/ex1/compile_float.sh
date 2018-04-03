#!/bin/bash

# The original source code folder of the mdatom code must be given as an argument
source_dir=$1

# Copy source dir here to be able to modify it freely
cp -r $source_dir src_float

# Replace all "double" by "float"
sed -i '' -e 's/double/float/g' $(find src_float/ -type f)

# Compile the project
mkdir build_float
cd build_float
cmake ../src_float
make
mv mdatom ../mdatom_float


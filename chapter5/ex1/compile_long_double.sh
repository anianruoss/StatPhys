#!/bin/bash

# The original source code folder of the mdatom code must be given as an argument
source_dir=$1

# Copy source dir here to be able to modify it freely
cp -r $source_dir src_long

# Replace all "double" by "long double"
sed -i '' -e 's/double/long double/g' $(find src_long/ -type f)

# Compile the project
mkdir build_long
cd build_long
cmake ../src_long
make
mv mdatom ../mdatom_long


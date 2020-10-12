#!/bin/bash

# Setup some safe shell options
set -eu -o pipefail
shopt -s nullglob

if [ ! -r run-tests.sh ]; then
    echo "$0: Must be run from the testing directory" >&2
    exit 1
fi

# Find all example cases
for d in ../examples/*/
do

	# Get case path
	casePath=${d%/}

	# Get case name
	caseName=${casePath##*/}

	# Print header
	printf "\nStoring new $caseName RefData!\n\n"

	# Clean/create the case directory
	rm -rf $caseName
	mkdir $caseName

	# Copy params.h from examples into src folder
	cp $casePath/params.h ../inc/params.h

	# Modifiy the write out frequencies for testing
	sed -i.bak "s|const int nSteps.*|const int nSteps = 500;|" ../inc/params.h
	sed -i.bak "s|const int tinfo.*|const int tinfo = nSteps / 50;|" ../inc/params.h
	sed -i.bak "s|const int tVTK.*|const int tVTK = nSteps / 10;|" ../inc/params.h
	sed -i.bak "s|const int tRestart.*|const int tRestart = nSteps / 5;|" ../inc/params.h

	rm -f ../inc/params.h.bak

	# Build LIFE
	(cd .. && make clean && make -j 8)

	# Create ref directory
	mkdir $caseName/RefData

	# Copy case to RefData
	cp ../LIFE $caseName/RefData/.
	cp ../inc/params.h $caseName/RefData/.

	# Check if there a geometry.config file
	if [ -d $casePath/input ]; then
		cp -r $casePath/input $caseName/RefData/.
	fi

	# Run the case
	(cd $caseName/RefData && ./LIFE)

	# If TurekHron case then run again to test restart feature
	if [ $caseName == "TurekHron" ]; then
		(cd $caseName/RefData && ./LIFE)
	fi

	# Print finish
	printf "Finished storing new $caseName RefData!\n\n"
done
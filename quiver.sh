#!/bin/bash

# This is a shell script that will create the input files for QUIVER,
# run it, and then calculate the KIEs based on the result.  It also
# cleans up the directory afterwards.
#
# usage:
# ./quiver.sh gs.out ts.out

# check command line arguments
if [ "$#" -ne 2 ]; then
    echo "usage: ./quiver.sh gs.out ts.out"
    exit 1
fi
if [ ! -f $1 ]; then
    echo "Could not find ground state file "$1"!"
    exit 1
fi
if [ ! -f $2 ]; then
    echo "Could not find transition state file "$2"!"
    exit 1
fi

# pre cleanup
rm -f *.q1 *.q2 *.q3 *.qout \ *

# prepare QUIVER input files
awk -f quiver_prep.awk quiver.config $1 $2
if [ $? -ne 0 ]; then
    echo "Something went wrong when preparing the input files."
    exit 1
fi

# run QUIVER on ground state
cp gs.q1 temp.q1
cp gs.q2 temp.q2
rm -f temp.qout
./quiver.exe &> gs.console
if [ $? -ne 0 ]; then
    echo "Something went wrong when running QUIVER on the ground state."
    exit 1
fi
cp temp.qout gs.qout

# run QUIVER on transition state
cp ts.q1 temp.q1
cp ts.q2 temp.q2
rm -f temp.qout
./quiver.exe &> ts.console
if [ $? -ne 0 ]; then
    echo "Something went wrong when running QUIVER on the ground state."
    exit 1
fi
cp temp.qout ts.qout

# run the analysis
awk -f quiver_analysis.awk gs.q3 gs.qout ts.qout
if [ $? -ne 0 ]; then
    echo "Something went wrong when performing the analysis."
    exit 1
fi

# clean up
rm -f *.q1 *.q2 *.q3 *.qout \ * *.console

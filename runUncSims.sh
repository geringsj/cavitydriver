#!/bin/bash

# for distributed smulation: run this script with the input config files as arguments.
# this script will start up to 4 single-threaded cavitydrivers with the 
# config-files from the arguments list, until no config files are left
while test $# -gt 0
do
	#echo $1
	if [ $(ps -e | grep cavitydriver | wc -l) -lt 4 ]
	then
		#echo "starting with $1"
		./../bin/cavitydriver $1 > /dev/null &
		shift
	fi
done


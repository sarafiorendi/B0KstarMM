#!/bin/sh

# First parameter: number of jobs
# Second parameter: q2 bin index

dir2Move="FLAFB10"
dir2Go="../tmp"

echo "Directory to move:" $dir2Move
echo "Directory to go:" $dir2Go

if [ $# -eq 0 ]
    then
    echo "Syapsis: source runMoveDirs.sh number_of_jobs q2bin"
else
    echo "Parameters:" $1 "and" $2
    
    i=1;
    while [ "$i" -le "$1" ]
      do
      echo "mv" $dir2Move"_"$i"_$2" $dir2Go
      mv $dir2Move"_"$i"_$2" $dir2Go
      i=$((i+1))
    done
fi

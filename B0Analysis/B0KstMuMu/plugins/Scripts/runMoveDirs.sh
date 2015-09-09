#!/bin/sh
# First parameter: number of jobs
# Second parameter: q2 bin index

echo "Parameters: " $1 " and " $2

i=1;
while [ "$i" -le "$1" ]
do
  echo "FLAFB05_"$i"_$2 tmp"
#  mv "FLAFB05_"$i"_$2" tmp
  i=$((i+1))
done

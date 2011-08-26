#!/bin/bash
# C = A - B
# bash colsubtract.sh A B > C
awk 'FNR==NR {arr[NR]=$1; next} { print $1-arr[FNR]}' $2 $1

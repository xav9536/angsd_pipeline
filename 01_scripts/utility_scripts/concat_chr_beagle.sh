#!/bin/bash

### This script will concatenate all per-chromosome beagle files from a folder into a single beagle file.

PREFIX=$1 ### This should include the folder containing the beagle files and the file name exclude "_chrXX"
SUFFIX="beagle" ### Can be changed to concatenate files other than beagle

## Create header
zcat $(ls -1 "$PREFIX"_chr[0-9][0-9]."$SUFFIX".gz | head -n1) | head -n1 | gzip > "$PREFIX".beagle.gz

## Concatenate file while skipping header
for f in $(ls -1 "$PREFIX"_chr[0-9][0-9]."$SUFFIX".gz)
    do
    zcat $f | tail -n +2 | gzip >> "$PREFIX"."$SUFFIX".gz
done

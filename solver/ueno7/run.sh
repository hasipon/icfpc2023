#!/bin/bash

for f in $(ls ../../problems.kyopro/*.kyopro);
do
    cat $f | ./a.out > $(basename $f .kyopro).stat 2>$(basename $f .kyopro).err &
done
wait

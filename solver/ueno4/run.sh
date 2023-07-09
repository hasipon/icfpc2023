#!/bin/bash

for f in $(ls ../../problems.kyopro/*.kyopro);
do
    cat $f | ./a.out | tee $(basename $f .kyopro).csv
done
wait

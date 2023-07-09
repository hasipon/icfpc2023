#!/bin/bash

for f in $(ls ../../problems.kyopro/*.kyopro);
do
    cat $f | ./a.out > ../../solutions/$(basename $f .kyopro)-ueno5.json 2>$(basename $f .kyopro).err &
done
wait

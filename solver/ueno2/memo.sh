#!/bin/bash

for i in `seq 1 55`
do
    echo $i.err:$(cat $i.err)
done

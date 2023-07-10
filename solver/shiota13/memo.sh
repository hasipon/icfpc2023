#!/bin/bash

for i in `seq 1 90`
do
    echo $i.err:$(cat $i.err)
done

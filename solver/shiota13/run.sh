#!/bin/bash

set -eux
cd $REPO_ROOT/solver/shiota13
g++ -Ofast -std=c++14 main.cpp
./a.out < $REPO_ROOT/problems.kyopro/${PROBLEM_ID}.kyopro > $REPO_ROOT/solutions/${PROBLEM_ID}-shiota13.json
cd shohei15-2 && cargo run --release

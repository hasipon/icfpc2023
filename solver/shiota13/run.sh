#! /bin/bash

set -eux
g++ -Ofast -std=c++14 main.cpp
./a.out > $REPO_ROOT/solutions/${PROBLEM_ID}-shiota13.json
cd shohei15-2 && cargo run --release

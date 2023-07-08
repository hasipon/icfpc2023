#!/bin/bash
set -x
SCRIPT_DIR=$(cd $(dirname $0) && pwd)

set -e
mkdir -p build
pushd build
  cmake .. && make
popd

set +e
for i in $(seq 1 25); do
PROBLEM_ID=${i} PROBLEM_PAM=../../../problems.pam/${PROBLEM_ID}.pam INPUT_ISL=../../../solutions.best/${PROBLEM_ID}.isl ./build/inada1
done

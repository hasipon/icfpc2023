#!/bin/bash

for i in $(seq 1 45); do \
  echo ${i}
  curl "https://api.icfpcontest.com/problem?problem_id=${i}" | jq -r .Success > "problem.json/${i}.json"
done

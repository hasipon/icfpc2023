set -eux
g++ -Ofast -std=c++14 main.cpp
for i in `seq 2 90`
do
	./a.out < ../../problems.kyopro/$i.kyopro > ../../solutions/${i}-shiota13.json
	PROBLEM_ID=$i REPO_ROOT=/home/shiota/icfpc2023 cargo run --release --manifest-path ./shohei15-2/Cargo.toml > ../../solutions/$i-shiota-13-v2.json &
done

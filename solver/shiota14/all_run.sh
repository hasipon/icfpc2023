set -eux
g++ -Ofast -std=c++14 main.cpp
cargo build --release --manifest-path ./shohei15-2/Cargo.toml
for i in `seq 1 1`
do
	./a.out < ../../problems.kyopro/$i.kyopro > ../../solutions/${i}-shiota14.tmp
	PROBLEM_ID=$i REPO_ROOT=/home/shiota/icfpc2023 ./shohei15-2/target/release/solver > ../../solutions/$i-shiota-14-v2.json &
done

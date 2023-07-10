for i in `seq 1 20`
do
	PROBLEM_ID=$i REPO_ROOT=/home/shiota/icfpc2023 /home/shiota/icfpc2023/solver/shiota13/run.sh > ../../solutions/$i-shiota-13-v2.json &
done

for i in `seq 1 90`
do
	./a.out < ../../problems.kyopro/$i.kyopro > $i-johniel10-2.json 2> $i-johniel10.err &
done

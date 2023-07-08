for i in `seq 1 55`
do
	./a.out < ../../problems.kyopro/$i.kyopro > $i-johniel6.json 2> $i-johniel6.err &
done

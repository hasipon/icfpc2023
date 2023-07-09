for i in `seq 1 90`
do
	./a.out < ../../problems.kyopro/$i.kyopro > $i-johniel9-2.json 2> $i-johniel9.err &
done

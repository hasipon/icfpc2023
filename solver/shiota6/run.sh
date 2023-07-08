for i in `seq 1 55`
do
	./a.out < ../../problems.kyopro/$i.kyopro > $i.out 2> $i.err &
done


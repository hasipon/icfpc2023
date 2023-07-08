# python3 conv.py ../../solutions/$1.json > $1.in
# ./a $1.in < ../../problems.kyopro/${1%%-*}.kyopro > $1.out 2> $1.err


for i in `seq 56 90`
do
	./a.out < ../../problems.kyopro/$i.kyopro > $i-johniel7.json 2> $i-johniel6.err &
done

python3 conv.py ../../solutions/$1.json > $1.in
./a $1.in < ../../problems.kyopro/${1%%-*}.kyopro > $1.out

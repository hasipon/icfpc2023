IDX=$(python3 get_b.py b$1.txt c$1.txt)
./c $1.out c$1.txt < ../../problems.kyopro/$1.kyopro > ../../solutions/$1-hasi10.$IDX.json

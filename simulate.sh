#!/bin/bash

P=( 0.45 0.455 0.46 0.465 0.47 0.475 0.48 0.485 0.49 )
k=$1
n=$2

for p in "${P[@]}"
do
	if [ ! -d "./data/p"$p"n"$2"" ]; then
		mkdir ./data/p"$p"n"$2"
	else
		continue
	fi

	for i in $(seq "$k")
	do
		./leath "$n" "$p" "$i" > ./data/p"$p"n"$n"/run"$i"p"$p"n"$n".data
		echo -ne 'Completed '"$i"' out of '"$k"' runs with probability: '"$p"'\r'
	done
done

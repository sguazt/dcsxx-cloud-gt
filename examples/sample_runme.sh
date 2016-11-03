#!/bin/sh

if [ "$#" -lt 2 ]; then
	echo "Usage: $0 <input scenario file> <output csv file>"
	exit 1
fi

scenario=$1
csv=$2

../src/sim	--csv $csv \
			--formation nash \
			--payoff shapley \
			--opt-relgap 0 \
			--scenario $scenario

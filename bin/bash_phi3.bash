#!/bin/bash
flies=13 # 13
conditions=2
taus=(4) #(4 8 16)
sets=1365
trials=8
for (( fly=4; fly<=$flies; fly++ )); do
	for (( condition=1; condition<=$conditions; condition++ )); do
		for tau in "${taus[@]}"; do
			for (( set=1; set<=$sets; set++ )); do
				for (( trial=1; trial<=$trials; trial++ )); do
					sbatch --job-name="f${fly}c${condition}tau${tau}s${set}t${trial}" --output="logs2/f${fly}c${condition}tau${tau}s${set}t${trial}.err" --error="logs/f${fly}c${condition}tau${tau}s${set}t${trial}.err" bash_phi3_sbatch.bash $fly $condition $tau $set $trial
					#python3 -i phi3.py $fly $set $condition $trial $tau
					#echo ./test.bash $fly $set $condition $trial $tau
				done
			done
		done
	done
done

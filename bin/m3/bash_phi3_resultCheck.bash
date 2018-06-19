#!/bin/bash
flies=10 # 13
conditions=2
taus=(8 16) #(4 8 16)
sets=1364
trials=8
prefix="split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_phithree_allPartitions_"
suffix=".mat"
for (( fly=1; fly<=$flies; fly++ )); do
	for (( condition=1; condition<=$conditions; condition++ )); do
		for tau in "${taus[@]}"; do
			for (( set=0; set<=$sets; set++ )); do
				for (( trial=1; trial<=$trials; trial++ )); do
					
					printf -v fly_padded "%02d" $fly
					printf -v set_padded "%04d" $set
					
					# Expected result file
					id="f${fly_padded}c${condition}tau${tau}s${set_padded}t${trial}"
					#results_file="results_split/fly${fly_padded}/${prefix}${id}${suffix}"
					results_file="results_split/${prefix}${id}${suffix}"
					
					# If result file doesn't exist, recompute
					if [ ! -e $results_file ]; then
						sbatch --job-name="${id}" --output="logs3/${id}.out" --error="logs3/${id}.err" bash_phi3_sbatch.bash $fly $condition $tau $(($set + 1)) $trial
						echo "${results_file}"
					fi
				done
			done
		done
	done
done


#!/bin/bash
flies=13 # 13
taus=(1) #(4 8 16)
nChannels=5
trials=8
channel_means_size=2
global_tpm=0
for (( fly=1; fly<=$flies; fly++ )); do
	for tau in "${taus[@]}"; do
		for (( trial=1; trial<=$trials; trial++ )); do
			squeue -u aleu6 > job_list
			jobs=$(wc -l < job_list)
			echo $jobs
			while [ $jobs -ge 16 ]; do
				sleep 30s
				squeue -u aleu6 > job_list
				jobs=$(wc -l < job_list)
				echo $jobs
				done
			sbatch --job-name="f${fly}tau${tau}t${trial}nCh${nChannels}g${global_tpm}" --output="logs/f${fly}tau${tau}t${trial}nCh${nChannels}g${global_tpm}.out" --error="logs/f${fly}tau${tau}t${trial}nCh${nChannels}g${global_tpm}.err" bash_phi3_chain_sbatch.bash $fly $tau $nChannels $trial $channel_means_size $global_tpm
		done
	done
done

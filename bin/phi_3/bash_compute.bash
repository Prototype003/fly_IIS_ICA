#!/bin/bash
nChannels=4
flies=1 # 13
conditions=(1 2)
set_ids=(1036)
taus=(1 2 3 4 8 12 16 24 32 48 64 128 256 512 4500 9000) #(4 8 16)
trials=1
global_tpm=1
tau_bin=1
for (( fly=1; fly<=$flies; fly++ )); do
	for condition in "${conditions[@]}"; do
		for set in "${set_ids[@]}"; do
			for tau in "${taus[@]}"; do
				for (( trial=1; trial<=$trials; trial++ )); do
					#for (( start_sample=1; start_sample<=$tau; start_sample++ )); do
					for (( start_sample=1; start_sample<=4; start_sample++ )); do
						squeue -u aleu6 > job_list
						jobs=$(wc -l < job_list)
						echo $jobs
						while [ $jobs -ge 35 ]; do
							sleep 30s
							squeue -u aleu6 > job_list
							jobs=$(wc -l < job_list)
							echo $jobs
						done
						sbatch --job-name="f${fly}tau${tau}t${trial}nCh${nChannels}g${global_tpm}" --output="logs/f${fly}tauBin${tau_bin}tau${tau}offset${start_sample}t${trial}nCh${nChannels}g${global_tpm}.out" --error="logs/f${fly}tauBin${tau_bin}tau${tau}offset${start_sample}t${trial}nCh${nChannels}g${global_tpm}.err" bash_compute_sbatch.bash $nChannels $fly $condition $set $tau $trial $global_tpm $tau_bin $start_sample
					done
				done
			done
		done
	done
done

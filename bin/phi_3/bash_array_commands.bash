#!/bin/bash
nChannels=4
flies=13 # 13
conditions=(1 2)
#set_ids=($(seq 1 1 1365)) # ($(seq first step last)) # (1036)
set_ids=($(seq 1 1 1))
taus=(4) # (1 2 3 4 8 12 16 24 32 48 64 128 256 512 4500 9000) #(4 8 16)
trials=8
global_tpm=0
tau_bin=0
sample_offset=0
> array_commands
for (( fly=1; fly<=$flies; fly++ )); do
	for condition in "${conditions[@]}"; do
		for set in "${set_ids[@]}"; do
			for tau in "${taus[@]}"; do
				for (( trial=1; trial<=$trials; trial++ )); do
					
					echo "python3 phi_compute.py $nChannels $fly $condition $set $tau $trial $global_tpm $tau_bin $sample_offset" >> array_commands
					
				done
			done
		done
	done
done

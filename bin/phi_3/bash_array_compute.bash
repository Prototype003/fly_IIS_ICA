#!/bin/bash
line_increment=20 # size of job array

# Fine number of parameter lines
lines=$(wc -l < array_commands) # Total number of jobs which need to be computed

# Loop through parameter lines
for (( line=1; line<=$lines; line=$line+$line_increment )); do
	squeue -u aleu6 > job_list
	jobs=$(wc -l < job_list)
	echo $jobs
	while [ $jobs -ge 11 ]; do
		sleep 60s
		squeue -u aleu6 > job_list
		jobs=$(wc -l < job_list)
		echo $jobs
	done
	sbatch --job-name="phi3_array" --output="logs/array_%A_%a.out" --error="logs/array_%A_%a.err" --array=1-$line_increment bash_array_sbatch.bash $line
done

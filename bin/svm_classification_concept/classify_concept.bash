#!/bin/bash

declare -a types=("part" "unpart" "both")

for type in "${types[@]}"; do
	command="matlab -nodisplay -nosplash -nodesktop -r \"main_classify_concept(${type}); exit\""
	sbatch --job-name="${type}" --output="logs/${type}.out" --error="logs/${type}.err" classify_concept_sbatch.bash "${command}"
done
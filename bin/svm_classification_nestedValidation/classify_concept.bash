#!/bin/bash

declare -a types=("part" "unpart" "both")
declare -a styles=("within" "across")

for style in "${styles[@]}"; do
	for type in "${types[@]}"; do
		command="matlab -nodisplay -nosplash -nodesktop -r \"main_classify_concept('${type}', '${style}'); exit\""
		sbatch --job-name="${style}${type}" --output="logs/${style}${type}.out" --error="logs/${style}${type}.err" classify_concept_sbatch.bash "${command}"
	done
done

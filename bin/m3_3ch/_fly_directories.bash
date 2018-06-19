#!/bin/bash
flies=13
for (( fly=1; fly<=$flies; fly++ )); do
	printf -v fly_padded "%02d" $fly
	mkdir "fly${fly_padded}"
	mv *_f"$fly_padded"*.mat fly"$fly_padded"/
done

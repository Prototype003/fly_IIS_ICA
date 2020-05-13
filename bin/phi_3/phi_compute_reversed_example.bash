#!/bin/bash

fly=1
set=876
tau=4

conds=2
trials=8

for (( cond=1; cond<=${conds}; cond=${cond}+1 )); do
	for (( trial=2; trial<=${trials}; trial=${trial}+1 )); do
		python phi_compute_reversed.py 4 ${fly} ${cond} ${set} ${tau} ${trial} 0 0 0
	done
done

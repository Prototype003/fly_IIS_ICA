#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Kai Xi,  Oct 2014
#              help@massive.org.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# $1: fly
# $2: condition
# $3: tau
# $4: nChannels
# $5: trial
# $6: channel_means_size
# $7: global_tpm
# Need to use variables OUTSIDE of this script, #SBATCH doesn't support variables: https://help.rc.ufl.edu/doc/Using_Variables_in_SLURM_Jobs
# SBATCH --job-name=f$1s$2c$3t$4t$5


# To set a project account for credit charging, 
#SBATCH --account=NCIfb7


# Request CPU resource for a serial job
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
# SBATCH --exclusive
#SBATCH --cpus-per-task=1

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-01:15:00


# To receive an email when job completes or fails
# SBATCH --mail-user=aleu6@student.monash.edu
# SBATCH --mail-type=END
# SBATCH --mail-type=FAIL


# Set the file for output (stdout)
# SBATCH --output=logs/myJob-fly2-%j.out

# Set the file for error log (stderr)
# SBATCH --error=logs/$1_$2_$3.err


# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name


# Command to run a serial job
module load python/3.6.0

python3 phi3_chain_tauBin.py $1 $2 $3 $4 $5 $6


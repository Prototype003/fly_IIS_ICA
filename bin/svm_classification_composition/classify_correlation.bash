#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Kai Xi,  Oct 2014
#              help@massive.org.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# $1: line counter
# Need to use variables OUTSIDE of this script, #SBATCH doesn't support variables: https://help.rc.ufl.edu/doc/Using_Variables_in_SLURM_Jobs
#SBATCH --job-name=correlation


# To set a project account for credit charging, 
#SBATCH --account=qb48


# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
# SBATCH --exclusive
#SBATCH --cpus-per-task=12

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-3:00:00

# SBATCH --qos=shortq

# To receive an email when job completes or fails
#SBATCH --mail-user=aleu6@student.monash.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the file for output (stdout)
#SBATCH --output=classify_coherence.out

# Set the file for error log (stderr)
#SBATCH --error=classify_coherence.err


# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name


# Job script

module load matlab/r2018a

#matlab -nodisplay -nosplash -nodesktop -r "2+2;exit"
matlab -nodisplay -nosplash -nodesktop -r "main_classify_correlation; exit"

# Installing pyphi should also install its dependencies
# Additional (not installed with pyphi) required packages: sklearn

import sys
import numpy as np
import scipy.signal as sp_signal
import scipy.stats as sp_stats
import pyphi
import pyphi_mimic
import itertools
from fly_phi import *

print("Packages imported")

# Setup ############################################################################

prep_detrend = 0;
prep_zscore = 0;

# Possible methods of binarisation (binarisation occurs all at once as preprocessing)
# binarise_global_median - discretise based on global median
# binarise_epoch_median - discretise based on local median within epoch
# binarise_global_pattern - discretise based on fitting patterns
# binarise_global_gradient - similar to pattern, except we generalise the shape to a general direction
binarise_method = "binarise_global_median"

# MATLAB (1:15) = PYTHON np.arange(1, 16)
# First parameter is inclusive
# Second parameter is exclusive (e.g. (x, y) gives from x to y-1)
# Third parameter of np.arange give sthe step size
flies = np.arange(0, 13)
nChannels = np.arange(2, 5)#(2, 16)
nBins = np.array([1]) # Script currently only supports the case of 1 bin
taus = np.array([4, 8, 16])
possible_partitions = np.array([None, None, 2, 6, 14])

nFlies = np.size(flies);

data_directory = "workspace_results/"
data_file = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim"

results_directory = "results/"
results_file = data_file + \
	"_detrend" + str(prep_detrend) + \
	"_zscore" + str(prep_zscore) + \
	"_nChannels" + str(nChannels[0]) + "t" + str(nChannels[-1]) + \
	"_tpmsSBN"

# Load data ############################################################################

loaded_data = load_mat(data_directory + data_file + ".mat")
fly_data = loaded_data['fly_data']
print("Fly data loaded")

#########################################################################################
# Remember that the data is in the form a matrix
# Matrix dimensions: sample(2250) x channel(15) x trial(8) x fly(13) x condition(2)

if prep_detrend == 1:
	fly_data = sp_signal.detrend(fly_data, 0)
if prep_zscore == 1:
	fly_data = sp_stats.zscore(fly_data, 0)

# Discretise
fly_data_discretised, n_values, channel_medians, *unused = globals()[binarise_method](fly_data)

# Get all channel combinations
channel_combinations = [] # List of lists (list elements are lists of differing lengths)
for channels in nChannels:
	channel_combinations.append(list(itertools.combinations(np.arange(fly_data.shape[1]), channels)))

# Results structure
tpms_sbn = [dict() for n in range(nChannels.size)]

# Following is repeated for each fly, bin count, tau lag, channel set, condition and trial (like in computation of phi-star)
for nChannels_counter in range(0, len(nChannels)):
	channel_sets = channel_combinations[nChannels_counter]
	partitions_n = possible_partitions[nChannels[nChannels_counter]]
	
	# Determine number of system states
	n_states = n_values ** len(channel_sets[0])
	
	# Store the number of channels and the associated channel sets, and the taus and nBins values
	tpms_sbn[nChannels_counter]['nChannels'] = nChannels[nChannels_counter]
	tpms_sbn[nChannels_counter]['channel_sets'] = channel_sets
	tpms_sbn[nChannels_counter]['taus'] = taus
	tpms_sbn[nChannels_counter]['nBins'] = nBins
	
	# Storage of TPMs (states x channels x channel sets x flies x conditions x taus)
	tpms = np.zeros((n_states, nChannels[nChannels_counter], len(channel_sets), fly_data.shape[3], fly_data.shape[4], len(taus)))
	
	for fly in flies:
		for condition in range(0, fly_data.shape[4]):
			for tau_counter in range(0, len(taus)):
				tau = taus[tau_counter]
				print ('fly' + str(fly) + ' condition' + str(condition) + ' tau' + str(tau) + ' channels_used' + str(len(channel_sets[0])))
				for channel_set_counter in range(0, len(channel_sets)):
					print('Channel set ' + str(channel_set_counter))
					sys.stdout.flush()
					channel_set = list(channel_sets[channel_set_counter]) # for indexing purposes
					
					# Build TPM for the set of channels
					# This is built across all trials
					# If it is to be built per trial, move this line to the 'for trial' loop
					tmp_data = fly_data_discretised[:, :, :, fly, condition] # There's an issue with vector indexing (combined with splicing and single value indexing), so can't use [:, channel_set, :, fly, condition] - have to separate into 2 lines
					tpm = build_tpm_sbn(tmp_data[:, channel_set, :], tau, n_values)
					tpms[:, :, channel_set_counter, fly, condition, tau_counter] = tpm
					
	tpms_sbn[nChannels_counter]['tpms'] = tpms

# Save ###########################################################################

save_mat(results_directory+results_file, {'tpms_sbn': tpms_sbn})
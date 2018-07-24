# Installing pyphi should also install its dependencies
# Additional (not installed with pyphi) required packages: sklearn

import os
import sys
sys.path.append('../')
import numpy as np
import scipy.signal as sp_signal
import scipy.stats as sp_stats
import pyphi
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
nChannels = np.arange(4, 5)#(2, 16)
nBins = np.array([1]) # Script currently only supports the case of 1 bin
taus = np.array([4, 8, 16])
possible_partitions = np.array([None, None, 2, 6, 14])

nFlies = np.size(flies);

data_directory = "../workspace_results/"
data_file = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim"

results_directory = "results/"
results_file = data_file + \
	"_detrend" + str(prep_detrend) + \
	"_zscore" + str(prep_zscore) + \
	"_nChannels" + str(nChannels[0]) + "t" + str(nChannels[-1]) + \
	"_phithree_allPartitions"

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

# Specify channel sets
# {5 6 7 8} for fly 1, awake, tau=4ms
# Remember - python indexing from 0, so subtract 1 from desired channel indices
a = 4
b = 5
c = 6
d = 7
outputs = ['ab.txt', 'ac.txt', 'ad.txt', 'bc.txt', 'bd.txt', 'cd.txt', 'abcd.txt']
channel_sets = [[a, b], [a, c], [a, d], [b, c], [b, d], [c, d], [a, b, c, d]]

condition = 0
fly = 1
tau = 4
set_counter = 1
for set_counter in range(0, len(channel_sets)):
	print(outputs[set_counter])
	channel_set = channel_sets[set_counter]
	n_states = n_values ** len(channel_set)
	
	# Get data
	tmp_data = fly_data_discretised[:, :, :, fly, condition] # samples x channels x trials
	# Flatten trials
	channel_data = tmp_data[:, :, 0]
	for trial in range(1, tmp_data.shape[2]):
		channel_data = np.concatenate((channel_data, tmp_data[:, :, trial]), axis=0)
	
	# Build TPM
	tpm = build_tpm(np.expand_dims(channel_data[:, channel_set], axis=2), tau, n_values)
	f = open(outputs[set_counter], 'w')
	f.write('state by state tpm:\n')
	f.write(str(tpm))
	f.write('\n\n')
	tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm)
	f.write('state by node tpm:\n');
	f.write(str(tpm_formatted))
	f.write('\n\n')
	
	# Build network and subsystem
	network = pyphi.Network(tpm_formatted)
	
	# Calculate phi for each possible state
	for state_index in range(0, n_states):
		# Figure out loli state
		state = pyphi.convert.loli_index2state(state_index, len(channel_set))
		print(state)
		
		# The subsystem is the full network
		subsystem = pyphi.Subsystem(network, state, network.node_indices)
		
		# Compute phi value
		big_mip = pyphi.compute.big_mip(subsystem)
		
		# Print
		f.write('state ' + str(state) + '\n')
		f.write(str(big_mip))
		f.write('\n\n');
	f.close()

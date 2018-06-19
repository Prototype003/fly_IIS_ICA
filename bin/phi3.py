# Installing pyphi should also install its dependencies
# Additional (not installed with pyphi) required packages: sklearn

import sys
import numpy as np
#import scipy.signal as sp_signal
#import scipy.stats as sp_stats
import pyphi
import pyphi_mimic
#import itertools
from fly_phi import *

print("Packages imported")

# Setup ############################################################################

fly = int(sys.argv[1]) # Only for loading data
condition = int(sys.argv[2]) # Only for loading data 
tau = int(sys.argv[3]) # Actual value
set = int(sys.argv[4]) - 1 # Python indexing from 0
trial = int(sys.argv[5]) # Only for loading data

n_values = 2 # binarised data
nPartitions = 14 # ways to bipartition 4 channels

data_directory = "fly_data_split/"
data_file_prefix = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim"
data_file_suffix = "_f" + str(fly) + "c" + str(condition) + "t" + str(trial)
data_file = data_file_prefix + data_file_suffix

# Load data ############################################################################

loaded_data = load_mat(data_directory + data_file + ".mat")
fly_data = loaded_data['data']
channel_sets = loaded_data['channel_sets'];
nChannels = loaded_data['nChannels'][0][0];

results_directory = "results_split/"
results_file_suffix = "_f" + "{0:0>2}".format(fly) + "c" + str(condition) + "tau" + str(tau) + "s" + "{0:0>4}".format(set) + "t" + str(trial)
results_file = data_file_prefix + \
	"_nChannels" + str(nChannels) + \
	"_phithree_allPartitions" + results_file_suffix

print("Data loaded")

#########################################################################################
# Remember that the data is in the form a matrix
# Matrix dimensions: sample(2250) x channel(15)

# Get channel set
channel_set = channel_sets[set, :]

# Determine number of system states
n_states = n_values ** len(channel_set)

# Results structure
phi = dict()
phi['nChannels'] = nChannels
phi['channel_set'] = channel_set
phi['tau'] = tau

# Initialise matrices
mips = np.empty(n_states, dtype=tuple)
state_counter = np.zeros(n_states)
state_phis = np.zeros(n_states)
state_partitions = np.empty((n_states, nPartitions), dtype=tuple)
state_partitions_phis = np.zeros((n_states, nPartitions))

# Build TPM
tpm = build_tpm(fly_data[:, channel_set, None], tau, n_values)
print("TPM built")

# Build the network and subsystem
# We are assuming full connection
tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm)
network = pyphi.Network(tpm_formatted)
print("Network built")

# Calculate all possible phi values (number of phi values is limited by the number of possible states)
for state_index in range(0, n_states):
	#print('State ' + str(state_index))
	# Figure out the state
	state = pyphi.convert.loli_index2state(state_index, nChannels)
	
	# As the network is already limited to the channel set, the subsystem would have the same nodes as the full network
	subsystem = pyphi.Subsystem(network, state, network.node_indices)
	
	#sys.exit()
	
	# Compute phi values for all partitions
	big_mip = pyphi_mimic.big_mip(subsystem)
	
	# Store phi and associated MIP
	state_phis[state_index] = big_mip[0].phi
	mips[state_index] = big_mip[0].cut
	
	# Store phis for all partition schemes
	for partition_counter in range(0, nPartitions):
		state_partitions_phis[state_index, partition_counter] = big_mip[1][partition_counter].phi
		state_partitions[state_index, partition_counter] = big_mip[1][partition_counter].cut
	
	print('State ' + str(state_index) + ' Phi=' + str(big_mip[0].phi))

# We average across the phi values previously computed for the channel set
# When averaging across samples, we weight by the number of times each state occurred
for sample_counter in range(0, fly_data.shape[0]):
	sample = fly_data[sample_counter, channel_set]
	
	# Determine the state
	state = pyphi.convert.state2loli_index(tuple(sample))
	
	# Add to state_counter
	state_counter[state] += 1
	
	# Add to state trajectory
	#state_trajectories[sample_counter, channel_set_counter, trial, fly, condition, tau_counter, bins_counter] = state
	
phi_total = 0
for state_index in range(0, len(state_counter)):
	if state_counter[state_index] > 0:
		# Add phi to total, weighted by the number of times the state occurred
		phi_total += state_phis[state_index] * state_counter[state_index]
phi_value = phi_total / np.sum(state_counter)


phi['phi'] = phi_value
phi['state_counters'] = state_counter
phi['state_phis'] = state_phis
phi['tpm'] = tpm
phi['mips'] = mips
phi['state_partitions'] = state_partitions
phi['state_partitions_phis'] = state_partitions_phis

# Save ###########################################################################

save_mat(results_directory+results_file, {'phi': phi})
print('saved ' + results_directory + results_file)
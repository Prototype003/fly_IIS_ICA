

import sys
sys.path.append('../')
import numpy as np
import scipy.signal as sp_signal
import scipy.stats as sp_stats
import pyphi
import itertools
from fly_phi import *

# Calculate for only 1 fly, 1 set size, loop using bash script

# Setup ############################################################################

# Parameters for loading data and calculating
fly = int(sys.argv[1]) # 1-13, only for loading and saving
tau = int(sys.argv[2]) # Actual tau in terms of how many samples to lag across (in ms, for 1000Hz sampling rate)
nChannels = int(sys.argv[3]) # Number of channels across which to compute phi
trial = int(sys.argv[4]) # 1-8, only for loading and saving
channel_mean_nChannels = int(sys.argv[5]) # Set size used to calculate channel mean values
channel_mean_nChannels = channel_mean_nChannels - 2 # 2 channels is in the first position (index 0)
global_tpm = int(sys.argv[6])

nConditions = 2 # awake/anest
n_values = 2 # binarised data

# Fly data location
data_directory = "../fly_data_split/"
data_file_prefix = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim"

# Channel means location
channel_data_directory = "../workspace_results/"
channel_data_file = "channelMeanPhis_globalTPM" + str(global_tpm) + "_tau" + str(4) + ".mat" # We always (currently) select based off of tau=4ms results

# Output location
results_directory = "results_phiChain_split/"
results_file_suffix = "_chainedPhi_f" + "{0:0>2}".format(fly) + "tau" + str(tau) + "t" + str(trial) 
results_file = data_file_prefix + \
	"_nChannels" + str(nChannels) + \
	"_globalTPM" + str(global_tpm) + \
	results_file_suffix

# Load data ############################################################################

# Raw data is loaded in the loop

# Channel phi means
loaded_data = load_mat(channel_data_directory + channel_data_file)
#channel_means = loaded_data['channel_means'][channel_mean_nChannels, 0]['phis'][0, 0][:, trial-1, fly-1, 0] # We want values based on wakeful channel means
channel_means = np.mean(loaded_data['channel_means'][channel_mean_nChannels, 0]['phis'][0, 0][:, :, fly-1, 0], axis=1) # Average phi values across trials, we want values based on wakeful channel means

# Get the channel set which produces the most phi
channel_set = np.argsort(channel_means)[-nChannels:] # Last nChannels channels after sorting by ascending channel means

#########################################################################################
# Remember that the data is in the form a matrix
# Matrix dimensions: sample(2250) x channel(15)

# Determine number of system states
n_states = n_values ** len(channel_set)

# Results structure
phi = dict()
phi['nChannels'] = nChannels
phi['channel_set'] = channel_set
phi['tau'] = tau

# Initialise matrices
phi_values = np.zeros(nConditions)
mips = np.empty((n_states, nConditions), dtype=tuple)
tpms = np.zeros((n_states, n_states, nConditions))
state_counters = np.zeros((n_states, nConditions))
state_phis = np.zeros((n_states, nConditions))

for condition in range(0, nConditions):
	
	# Load data
	data_file_suffix = "_f" + str(fly) + "c" + str(condition+1) + "t" + str(trial)
	data_file = data_file_prefix + data_file_suffix + ".mat"
	loaded_data = load_mat(data_directory + data_file)
	fly_data = loaded_data['data']
	print("Data loaded")
	
	# Build TPM
	tpm = build_tpm(fly_data[:, channel_set, None], tau, n_values)
	tpms[:, :, condition] = tpm
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
		big_mip = pyphi.compute.big_mip(subsystem)
		
		# Store phi and associated MIP
		state_phis[state_index, condition] = big_mip.phi
		mips[state_index, condition] = big_mip.cut
		
		print('State ' + str(state_index) + ' Phi=' + str(big_mip.phi))

	# We average across the phi values previously computed for the channel set
	# When averaging across samples, we weight by the number of times each state occurred
	for sample_counter in range(0, fly_data.shape[0]):
		sample = fly_data[sample_counter, channel_set]
		
		# Determine the state
		state = pyphi.convert.state2loli_index(tuple(sample))
		
		# Add to state_counter
		state_counters[state, condition] += 1
		
		# Add to state trajectory
		#state_trajectories[sample_counter, channel_set_counter, trial, fly, condition, tau_counter, bins_counter] = state
		
	phi_total = 0
	for state_index in range(0, n_states):
		if state_counters[state_index, condition] > 0:
			# Add phi to total, weighted by the number of times the state occurred
			phi_total += state_phis[state_index, condition] * state_counters[state_index, condition]
	phi_values[condition] = phi_total / np.sum(state_counters[:, condition])


phi['phis'] = phi_values
phi['state_counters'] = state_counters
phi['state_phis'] = state_phis
phi['tpms'] = tpms
phi['mips'] = mips

# Save ###########################################################################

save_mat(results_directory+results_file, {'phi': phi})
print('saved ' + results_directory + results_file)
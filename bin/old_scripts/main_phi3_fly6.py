# Installing pyphi should also install its dependencies
# Additional (not installed with pyphi) required packages: sklearn

import sys
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
flies = np.array([5])
nChannels = np.arange(4, 5)#(2, 16)
nBins = np.array([1]) # Script currently only supports the case of 1 bin
taus = np.array([4, 8, 16])

nFlies = np.size(flies);

data_directory = "workspace_results/"
data_file = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim"

results_directory = "results/"
results_file = data_file + \
	"_detrend" + str(prep_detrend) + \
	"_zscore" + str(prep_zscore) + \
	"_nChannels" + str(nChannels[0]) + "t" + str(nChannels[-1]) + \
	"_phithree" + "fly" + str(flies[0])

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
phis = [dict() for n in range(nChannels.size)]

# Following is repeated for each fly, bin count, tau lag, channel set, condition and trial (like in computation of phi-star)
for nChannels_counter in range(0, len(nChannels)):
	channel_sets = channel_combinations[nChannels_counter]
	
	# Determine number of system states
	n_states = n_values ** len(channel_sets[0])
	
	# Store the number of channels and the associated channel sets, and the taus and nBins values
	phis[nChannels_counter]['nChannels'] = nChannels[nChannels_counter]
	phis[nChannels_counter]['channel_sets'] = channel_sets
	phis[nChannels_counter]['taus'] = taus
	phis[nChannels_counter]['nBins'] = nBins
	phis[nChannels_counter]['mips'] = []
	
	# Initialise phi_matrix (channel sets x trials x flies x conditions x taus x nBins)
	phi_threes = np.zeros((len(channel_sets), fly_data.shape[2], 1, fly_data.shape[4], len(taus), len(nBins)))
	
	# Storage of state counts (states x channel sets x trials x flies x conditions x taus x nBins) and phis (-trials and -nBins dimensions)
	state_counter = np.zeros((n_states, len(channel_sets), fly_data.shape[2], 1, fly_data.shape[4], len(taus), len(nBins)))
	state_phis = np.zeros((n_states, len(channel_sets), 1, fly_data.shape[4], len(taus), len(nBins))) # This will hold the phi values of each state
	
	# Storage of state trajectory (samples x channel sets x trials x flies x conditions x taus x nBins)
	#state_trajectories = np.zeros((fly_data.shape[0], len(channel_sets), fly_data.shape[2], 1, fly_data.shape[4], len(taus), len(nBins)))
	
	# Storage of TPMs (states x states x channel sets x flies x conditions x taus)
	tpms = np.zeros((n_states, n_states, len(channel_sets), 1, fly_data.shape[4], len(taus)))
	
	for fly_counter in range(0, len(flies)):
		fly = flies[fly_counter]
		for condition in range(0, fly_data.shape[4]):
			for tau_counter in range(0, len(taus)):
				tau = taus[tau_counter]
				print ('fly' + str(fly) + ' condition' + str(condition) + ' tau' + str(tau) + ' channels_used' + str(len(channel_sets[0])))
				for channel_set_counter in range(0, len(channel_sets)):
					print('Channel set ' + str(channel_set_counter))
					channel_set = list(channel_sets[channel_set_counter]) # for indexing purposes
					
					# Build TPM for the set of channels
					# This is built across all trials
					# If it is to be built per trial, move this line to the 'for trial' loop
					tmp_data = fly_data_discretised[:, :, :, fly, condition] # There's an issue with vector indexing (combined with splicing and single value indexing), so can't use [:, channel_set, :, fly, condition] - have to separate into 2 lines
					tpm = build_tpm(tmp_data[:, channel_set, :], tau, n_values)
					tpms[:, :, channel_set_counter, fly_counter, condition, tau_counter] = tpm
					
					# Build the network and subsystem
					# We are assuming full connection
					tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm)
					network = pyphi.Network(tpm_formatted)
					
					# Calculate all possible phi values (number of phi values is limited by the number of possible states)
					for state_index in range(0, n_states):
						print('State ' + str(state_index))
						# Figure out the state
						state = pyphi.convert.loli_index2state(state_index, nChannels[nChannels_counter])
						
						# As the network is already limited to the channel set, the subsystem would have the same nodes as the full network
						subsystem = pyphi.Subsystem(network, state, network.node_indices)
						
						# Compute phi (via mip)
						big_mip = pyphi.compute.big_mip(subsystem)
						state_phis[state_index, channel_set_counter, fly_counter, condition, tau_counter] = big_mip.phi
						phis[nChannels_counter]['mips'].append(big_mip.cut)
					
					for bins_counter in range(0, len(nBins)):
						bins = nBins[bins_counter]
						for trial in range(0, fly_data.shape[2]):
							
							# We average across the phi values previously computed for the channel set
							# When averaging across samples, we weight by the number of times each state occurred
							
							# Count the number of times each state occurs within the trial
							for sample_counter in range(0, fly_data.shape[0]):
								sample = fly_data_discretised[sample_counter, channel_set, trial, fly, condition]
								
								# Determine the state
								state = pyphi.convert.state2loli_index(tuple(sample))
								
								# Add to state_counter
								state_counter[state, channel_set_counter, trial, fly_counter, condition, tau_counter, bins_counter] += 1
								
								# Add to state trajectory
								#state_trajectories[sample_counter, channel_set_counter, trial, fly_counter, condition, tau_counter, bins_counter] = state
								
							phi_total = 0
							for state_index in range(0, len(state_counter)):
								if state_counter[state_index, channel_set_counter, trial, fly_counter, condition, tau_counter, bins_counter] > 0:
									# Add phi to total, weighted by the number of times the state occurred
									phi_total += state_phis[state_index, channel_set_counter, fly_counter, condition, tau_counter, bins_counter] * state_counter[state_index, channel_set_counter, trial, fly_counter, condition, tau_counter, bins_counter]
							phi_threes[channel_set_counter, trial, fly_counter, condition, tau_counter, bins_counter] = phi_total / np.sum(state_counter[:, channel_set_counter, trial, fly, condition, tau_counter, bins_counter])
	phis[nChannels_counter]['phi_threes'] = phi_threes
	phis[nChannels_counter]['state_counters'] = state_counter
	phis[nChannels_counter]['state_phis'] = state_phis
	#phis[nChannels_counter]['state_trajectories'] = state_trajectories
	phis[nChannels_counter]['tpms'] = tpms

# Save ###########################################################################

save_mat(results_directory+results_file, {'phis': phis})

# # Old code #######################################################################
# # This builds the entire TPM, then marginalises out irrelevant nodes
# # This takes too long because of network validation by pyphi
# # Thus it's better to build separate TPMs (and networks) for each channel set, where the subsystem is the whole network
# # rather than building the whole network and then limiting it to the relevant subsystem
# # Because of the long running time, sections after (and during) network creation have not been properly tested
# for fly in flies:
	# for condition in range(0, fly_data.shape[4]):
		# for tau in taus:
			# # Build transition probability matrix (we can build the whole matrix, for all channels, once, then later extract the relevant parts of it)
			# # For now we will build across all trials, for all channels
			# # We can marginalise out irrelevant channels (nodes) using pyphi.utils.marginalize_out
			# tpm_whole = build_tpm(fly_data_discretised[:, :, :, fly, condition], tau, n_values)
			# print('tpm built for fly' + str(fly) + ' condition ' + str(condition) + ' tau ' + str(tau))
			# channels_all = np.arange(0, fly_data.shape[1])
			# network = pyphi.Network(tpm_whole)
			# print('network built')
			# for nChannels_counter in range(0, len(nChannels)):
				# channel_sets = channel_combinations[nChannels_counter]
				# print('Channels in subsystem: ' + str(len(channel_sets[0])))
				# for channel_set in channel_sets:
					
					# channel_set = list(channel_set) # For indexing purposes
					# """ Unnecessary because this should be done by pyphi itself, when setting the subsystem
					# # Find channels which are not included in the current set
					# unincluded_channels = np.setdiff1d(channels_all, channel_set)
					
					# # Marginalise out unincluded channels from TPM
					# # Current version of pyphi only removes a single node at a time (future version will be able to remove a list of nodes)
					# # Seems like marginalize_out takes a state-by-node TPM: http://pyphi.readthedocs.io/en/stable/examples/conditional_independence.html?highlight=margin
					# tpm = pyphi.convert.state_by_state2state_by_node(np.copy(tpm_whole))
					# for channel in unincluded_channels:
						# tpm = pyphi.utils.marginalize_out(channel, tpm)
					# tpm = np.squeeze(tpm[:, :, :, channel_set]) # Restrict purview to relevant channels
					# """
					
					# # Determine number of system states
					# n_states = n_values ** len(channel_set)
					
					# for bins in nBins:
						# for trial in range(0, fly_data.shape[2]):
							
							# # To avoid calculating phi for every sample, we calculate it for each state which has occurred (e.g. 2^8 calculations for 8 channels, instead of 2250)
							# # Then when averaging across samples, we weight by the number of times each state has occurred
							
							# # Store how many times each state occurs within a trial
							# state_counter = np.zeros((n_states, 1))
							
							# for sample_counter in range(0, fly_data.shape[0]):
								# sample = fly_data[sample_counter, channel_set, trial, fly, condition]
								
								# # Determine the state
								# state = state2loli_index(tuple(sample))
								
								# # Add to state_counter
								# state_counter[state] += 1
								
							# # Calculate phi for each state which occurred
							# phis = np.zeros((n_states, 1))
							# phi_total = 0
							# for state_index in range(0, len(state_counter)):
								# if frequency > 0:
									# # Figure out the state
									# state = loli_index2state(state_index)
									# # Pad in other nodes (states given to Subsystem are of the network)
									# # I'm assuming that nodes not in channel_set are marginalised out, such that the state of these nodes is irrelevant
									# network_state = np.zeros((1, fly_data.shape[1]))
									# channel_counter = 0
									# for channel in channel_set:
										# network_state[channel] = state[channel_counter]
										# channel_counter += 1
									# # Limit the network to the set of channels being considered, and provide the relevant state
									# subsystem = pyphi.Subsystem(network, tuple(network_state), channel_set)
									# # Compute phi
									# phis[state] = pyphi.compute.big_phi(subsystem)
									# phi_total += phis[state] * state_counter[state_index]
							# phi = phi_total / np.sum(state_counter[:, channel_set_counter, trial, fly, condition, tau_counter, bins_counter])
							# print(phi)

import os
import sys
sys.path.append('../')
import numpy as np
import scipy.signal as sp_signal
import scipy.stats as sp_stats
import pyphi
import itertools
from fly_phi import *

# Calculate for only 1 fly, 1 set size, etc (see input parameters below) - loop using bash script
# Get all values - all mechanisms/purviews, and all partitions for each

# Modification of pyphi.subsystem.find_mip, to return values for all bipartitions
def find_mip(subsystem, direction, mechanism, purview):
	"""Return the minimum information partition for a mechanism over a
	purview.

	Args:
		direction (str): Either |past| or |future|.
		mechanism (tuple(int)): The nodes in the mechanism.
		purview (tuple(int)): The nodes in the purview.

	Returns:
		mip (|Mip|): The mininum-information partition in one temporal
			direction.
	"""
	repertoire = subsystem._get_repertoire(direction)

	# We default to the null MIP (the MIP of a reducible mechanism)
	mip = pyphi.subsystem._null_mip(direction, mechanism, purview)

	#if not purview:
	#	return mip

	phi_min = float('inf')
	# Calculate the unpartitioned repertoire to compare against the
	# partitioned ones
	unpartitioned_repertoire = repertoire(mechanism, purview)
	
	bipartitions = pyphi.subsystem.mip_bipartitions(mechanism, purview)
	
	partitions = [] # stores the results for a given partition
	bipartition_list = [] # stores the bipartition schemes (the parts)
	
	# Loop over possible MIP bipartitions
	for part0, part1 in bipartitions:
		# Find the distance between the unpartitioned repertoire and
		# the product of the repertoires of the two parts, e.g.
		#   D( p(ABC/ABC) || p(AC/C) * p(B/AB) )
		part1rep = repertoire(part0.mechanism, part0.purview)
		part2rep = repertoire(part1.mechanism, part1.purview)
		partitioned_repertoire = part1rep * part2rep

		phi = pyphi.utils.hamming_emd(unpartitioned_repertoire,
								partitioned_repertoire)
		phi = round(phi, pyphi.config.PRECISION)
		
		# Add bipartition results to list
		partition = pyphi.subsystem.Mip(direction=direction,
			mechanism=mechanism,
			purview=purview,
			partition=(part0, part1),
			unpartitioned_repertoire=unpartitioned_repertoire,
			partitioned_repertoire=partitioned_repertoire,
			phi=phi)
		
		partitions.append(partition)#.to_json())
		bipartition_list.append((part0, part1))

		# Update MIP if it's more minimal.
		if phi < phi_min:
			phi_min = phi
			# TODO: Use properties here to infer mechanism and purview from
			# partition yet access them with `.mechanism` and `.purview`.
			mip = partition
	print(partitions)
	print(bipartitions)
	print(mip)
	return partitions, bipartitions, mip

# Setup ############################################################################

# Parameters for loading data and calculating
nChannels = int(sys.argv[1]) # Number of channels across which to compute phi
fly = int(sys.argv[2]) # 1-indexed
condition = int(sys.argv[3]) # 1-indexed
set = int(sys.argv[4]) # 1-indexed
tau = int(sys.argv[5]) # Actual tau in terms of how many samples to lag across (in ms, for 1000Hz sampling rate)
trial = int(sys.argv[6]) # 1-indexed
global_tpm = int(sys.argv[7]) # 0=8 trials (2250 samples); 1=1 trial (18000 samples)
tau_bin = int(sys.argv[8]) # 0=don't average across tau samples, use stepsize tau; 1=average across tau samples, afterwards use stepsize 1
sample_offsets = int(sys.argv[9]) # Compute TPM across sample offsets before binning (for when tau_bin == 1)

# Fly data location
data_directory = "../workspace_results/"
data_file_prefix = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim"
data_file = data_file_prefix + ".mat"

# Output location
results_directory = "results_split/"
if not os.path.exists(results_directory):
	os.makedirs(results_directory)

# tau string for results file
if tau_bin == 1:
	tau_type = "tauBin"
	tau_string = tau_type + str(tau) + "binOffset" + str(sample_offsets)
else:
	tau_type = "tau"
	tau_string = tau_type + str(tau)

# Results file
results_file_suffix = "_nChannels" + str(nChannels) + "_globalTPM" + str(global_tpm) + "_f" + "{0:0>2}".format(fly) + "c" + str(condition) + tau_string + "s" + "{0:0>4}".format(set) + "t" + str(trial)
results_file = data_file_prefix + results_file_suffix + "_example.mat"

# Load data ############################################################################

loaded_data = load_mat(data_directory + data_file)
fly_data = loaded_data['fly_data']
print("Fly data loaded")

# Filter for data from channel set of specific fly
channel_combinations = list(itertools.combinations(np.arange(fly_data.shape[1]), nChannels)) # All channel combinations
channel_set = np.asarray(channel_combinations[set-1])

fly_data = fly_data[:, :, :, fly-1, condition-1] # From here, only holds data for specified fly, condition
fly_data = fly_data[:, channel_set, :] # Data for channel set (can't do with fly,condition because dimension order changes for some reason)

# Get trial
if global_tpm == 1:
	# Flatten trials if using global TPM
	flattened_data = fly_data[:, :, 0]
	for trial in range(1, fly_data.shape[2]):
		flattened_data = np.concatenate((flattened_data, fly_data[:, :, trial]), axis=0)
	fly_data = flattened_data
else:
	# Get specific trial
	fly_data = fly_data[:, :, trial-1] # From here, only holds data for specified trial

print("Specific data obtained")

# Preprocess ############################################################################

tpm_built = 0
if tau_bin == 1:
	if sample_offsets == 1:
		n_values = 2
		tpm, transition_counters = build_tpm_bin_offsets(fly_data, n_values, tau)
		tpm_formatted = tpm
		tpm_built = 1
	else:
		# Downsample by averaging in bins of length tau
		fly_data = tau_resample(fly_data[:, :, None, None, None], tau)
		fly_data = fly_data[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions
		tau_step = 1
else:
	tau_step = tau

if tpm_built == 0:
	# Binarise by median split
	fly_data, n_values, medians = binarise_trial_median(fly_data[:, :, None, None, None])
	fly_data = fly_data[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions

	# Build TPM
	tpm = build_tpm(fly_data[:, :, None], tau_step, n_values)
	tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm)

# Build sbn TPM for whole system and for split system
tpm_sbn = build_tpm_sbn(fly_data[:, :, None], tau_step, n_values)
tpm_sbn_0 = build_tpm_sbn(fly_data[:, 0, None, None], tau_step, n_values) # Need to inject extra dimension as the channels dimension is collapsed for one channel
tpm_sbn_1 = build_tpm_sbn(fly_data[:, 1, None, None], tau_step, n_values)
# Build split system TPM
tpm_sbn_ind = [
	[tpm_sbn_0[0][0], tpm_sbn_1[0][0]],
	[tpm_sbn_0[1][0], tpm_sbn_1[0][0]],
	[tpm_sbn_0[0][0], tpm_sbn_1[1][0]],
	[tpm_sbn_0[1][0], tpm_sbn_1[1][0]]] # Hardcoded for 2ch
	
print("TPM built")

# if sample_offsets == 1:
	# n_values = 2
	# tpm, transition_counters = build_tpm_bin_offsets(fly_data, n_values, tau)
	# tpm_formatted = tpm
# else:
	# if tau_bin == 1:
		# # Downsample by averaging in bins of length tau
		# fly_data = tau_resample(fly_data[:, :, None, None, None], tau)
		# fly_data = fly_data[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions
		# tau_step = 1
	# else:
		# tau_step = tau
	
	# # Binarise by median split
	# fly_data, n_values, medians = binarise_trial_median(fly_data[:, :, None, None, None])
	# fly_data = fly_data[:, :, 0, 0, 0] # Get rid of those appended singleton dimensions

	# # Build TPM
	# tpm = build_tpm(fly_data[:, :, None], tau_step, n_values)
	# tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm)
	# #tpm = build_tpm(np.flip(fly_data[:, :, None], 0), tau_step, n_values)
	# #tpm, tmp = build_tpm_sbn_normalise(fly_data[:, :, None], tau_step, n_values, 9000)
# print("TPM built")

# Build the network and subsystem
# We are assuming full connection
network = pyphi.Network(tpm_formatted)
print("Network built")

#########################################################################################
# Remember that the data is in the form a matrix
# Matrix dimensions: sample(2250) x channel(15)

# Determine number of system states
n_states = n_values ** len(channel_set)
n_states = 1 # we only want an example for one state

# Results structure
phi = dict()
phi['nChannels'] = nChannels
phi['channel_set'] = channel_set + 1 # Save 1-indexed values
phi['tau'] = tau

# Initialise results storage structures
phi_value = 0; # Doesn't need initialisation
mips = np.empty((n_states), dtype=tuple)
big_mips = np.empty((n_states), dtype=object)
state_counters = np.zeros((n_states))
state_phis = np.zeros((n_states))

# sys.exit()

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
	
	# Get all partition values for all mechanisms/purviews
	mech_list_unpartitioned = []
	mech_list_partitioned = []
	
	# Storage of the results for N nodes
	# number of mechanisms = powerset
	# number of purviews = powerset
	# number of mechanisms = number of purviews
	# vector of mechanisms
	# vector of purviews
	# matrix (mechanisms x purviews) of partitions
	# number of partitions depends on mechanism/purview
	
	# For each mechanism in the full system
	mechanisms = pyphi.utils.powerset(subsystem.node_indices)
	for mechanism in mechanisms:
		
		purviews = pyphi.utils.powerset(subsystem.node_indices)
		purview_phis = []
		
		# For each purview
		for purview in purviews:
			
			partition_phis, partitions, mip = find_mip(subsystem, 'past', mechanism, purview)
			partition_phis_future, partitions_future, mip_future = find_mip(subsystem, 'future', mechanism, purview)
			
			# Concatenate past and future
			partition_phis.extend(partition_phis_future)
			partitions.extend(partitions_future)
			mip_past_future = {'past' : mip, 'future' : mip_future}
			
			# Add to list
			purview_phis.append({'purview' : purview, 'partitions' : partitions, 'partition_phis' : partition_phis, 'mip' : mip_past_future})
			
		# Add mechanism and purview results to list
		mech_list_unpartitioned.append({'mechanism' : mechanism, 'purview_phis' : purview_phis})
	
	# For each mechanism in the partitioned system
	subsystem_partitioned = pyphi.Subsystem(network, state, network.node_indices, cut=big_mip.cut)
	mechanisms = pyphi.utils.powerset(subsystem.node_indices)
	for mechanism in mechanisms:
		
		purviews = pyphi.utils.powerset(subsystem_partitioned.node_indices)
		purview_phis = []
		
		# For each purview
		for purview in purviews:
			
			partition_phis, partitions, mip = find_mip(subsystem_partitioned, 'past', mechanism, purview)
			partition_phis_future, partitions_future, mip_future = find_mip(subsystem_partitioned, 'future', mechanism, purview)
			
			# Concatenate past and future
			partition_phis.extend(partition_phis_future)
			partitions.extend(partitions_future)
			mip_past_future = {'past' : mip, 'future' : mip_future}
			
			# Add to list
			purview_phis.append({'purview' : purview, 'partitions' : partitions, 'partition_phis' : partition_phis, 'mip' : mip_past_future})
			
		# Add mechanism and purview results to list
		mech_list_partitioned.append({'mechanism' : mechanism, 'purview_phis' : purview_phis})
	
	#sys.exit()
	
	# Store phi and associated MIP
	state_phis[state_index] = big_mip.phi
	mips[state_index] = big_mip.cut
	
	# MATLAB friendly storage format (python saves json as nested dict)
	big_mip_json = big_mip.to_json()
	# Sort each constellation by their mechanisms
	big_mip_json['partitioned_constellation'] = sort_constellation_by_mechanism(big_mip_json['partitioned_constellation'])
	big_mip_json['unpartitioned_constellation'] = sort_constellation_by_mechanism(big_mip_json['unpartitioned_constellation'])
	# Store big_mip
	big_mips[state_index] = big_mip_json
	
	print('State ' + str(state_index) + ' Phi=' + str(big_mip.phi))

# # We average across the phi values previously computed for the channel set
# # When averaging across samples, we weight by the number of times each state occurred
# if sample_offsets == 0:
	# for sample_counter in range(0, fly_data.shape[0]):
		# sample = fly_data[sample_counter, :]
		
		# # Determine the state
		# state = pyphi.convert.state2loli_index(tuple(sample))
		
		# # Add to state_counter
		# state_counters[state] += 1
# else:
	# state_counters = transition_counters

# phi_total = 0
# for state_index in range(0, n_states):
	# if state_counters[state_index] > 0:
		# # Add phi to total, weighted by the number of times the state occurred
		# phi_total += state_phis[state_index] * state_counters[state_index]
# phi_value = phi_total / np.sum(state_counters)

#phi['phi'] = phi_value
phi['state_counters'] = state_counters
phi['big_mips'] = big_mips
phi['state_phis'] = state_phis
phi['tpm'] = tpm
phi['mips'] = mips

# Save ###########################################################################
#sys.exit()
save_mat(results_file,
	{
		'phi': phi,
		'tpm_sbn': tpm_sbn,
		'tpm_sbn_0': tpm_sbn_0,
		'tpm_sbn_1': tpm_sbn_1,
		'tpm_sbn_ind': tpm_sbn_ind,
		'tpm_sbs': tpm,
		'big_mip': big_mip_json,
		'mech_list_unpartitioned': mech_list_unpartitioned,
		'mech_list_partitioned': mech_list_partitioned})
#save_mat(results_directory+results_file, {'phi': phi})

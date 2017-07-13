"""
To find package install locations, use:
import site; site.getsitepackages()

Anaconda (not standardalone Python):
Location should be "C:/Users/this_/Anaconda3/Lib/site-packages"
pyphi is installed in the above directory
"""

print("Setting up environment")

import pyphi as pyphi
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import sklearn.preprocessing as skp
import copy
import itertools

print("packages imported")




data_directory = "../../flies/fly_data/trials_anesthDescAdded11092014_bPlrRerefTyp2/Analyzed_WalkingDrosDror115610072014/" # 8 channels (after bipolar re-referencing - need confirmation, should create function for re-referencing)
#data_directory = "../../flies/fly_data/trials_anesthDescAdded11092014/Analyzed_WalkingDrosDror115610072014/" # 16 channels (raw)

"""
Data assumptions:

trials.mat holds a single Matlab variable - a struct array "trials"
trials contains the fields onTrial, startTime, endTime, LFP, duration, nFlicker, flickerRate, flickerRate2, anesthLabel, anesthNm



Accessing the fields via Matlab:
	foo = [trials.onTrial] (returns 1 x N vector)
	foo = [trials.LFP] (returns n_channels x n_samples 2D matrix; matrix is concatenation of matrices along the 1st (channels) dimension)
"""

loaded_fly_data = sio.loadmat(data_directory + "trials.mat")
fly_data = loaded_fly_data["trials"]
print("Data format: " + str(type(fly_data)) + " : " + str(fly_data.shape))

lfp = fly_data["LFP"][0, :] # LFP is a 1D vector (but numpy has it as a 2D matrix with the size of the first dimension=1) [0, :] gets the first (in this case, only) element in the 1st dimension, and all in the 2nd dimension
anesth_label = fly_data["anesthLabel"][0, :]
on_trial = fly_data["onTrial"][0, :]


"""
Time series must correspond to the length of a trial - when taking all trials, we likely can better populate the TPM (n_trials per channel, per timepoint)

If the time series is across the whole experiment, then the number of samples with which to build the TPM becomes much smaller (1 per channel, per timepoint)

Maybe can build TPM using all trials, then apply the same TPM to all samples / time windows
Is this valid if stimulus is constant?
What are the effects of 2.3s stim on / 0.8s stim off repetition?
"""


# Get all rests (times between on trials) for the air condition
lfp_air = lfp[np.ix_(np.logical_and(anesth_label==1, on_trial==0))]


def concat_arrays_in_array(matrix_array):
	"""
	Concatenates each matrix (numpy array, as stored in a numpy array) into a large matrix along the 1st dimension (number of elements in the 1st dimension stays the same; concatenates the other dimensions)
	Inputs:
		numpy array of numpy arrays - the first dimension of these arrays should be constant
			e.g. matrix_array[0] gives the first array
			e.g. matrix_array[1] gives the second array, etc.
	Outputs:
	"""
	
	# Take first trial
	concatenated = matrix_array[0]
	
	# Concatenate remaining trials
	for matrix in range(1, matrix_array.shape[0]):
		concatenated = np.concatenate((concatenated, matrix_array[matrix]), axis=1)
	
	return concatenated

# Concatenate all trials together (for finding the global median - for binarisation)
lfp_air_concatenated = concat_arrays_in_array(lfp_air)

"""
Binarisation:

All samples should be binarised
Median-split - get global median; if sample is greater than median, set to 1, else 0

"""

# Find global median per channel (using concatenated trials)
channel_medians = np.median(lfp_air_concatenated, axis=1) # Median across samples, per channel

def threshold_binarise_arrays_in_array(matrix_array, thresholds):
	"""
	Binarises matrices of an array, using the provided thresholds - each threshold corresponds to each element in the first dimension
	Each element in the second dimension is binarised using the threshold for their position in the first dimension
	Inputs:
		matrix_array - array of 2D matrices, all with the same number of elements in the first dimension
		thresholds - 1D vector which holds the thresholds to binarise by
			Should have the same number of elements as the first dimensions of the matrices
	Outputs:
	"""
	import sklearn.preprocessing as skp
	
	matrix_array = copy.deepcopy(matrix_array) # We don't want affect scope outside the function
	
	for matrix in range(matrix_array.shape[0]): # i.e. for each matrix
		for threshold_dimension in range(thresholds.shape[0]): # i.e. for each row (threshold_dimension is a D1 (row) counter)
			matrix_array[matrix][threshold_dimension, :] = skp.binarize(matrix_array[matrix][threshold_dimension, :].reshape(1, -1), thresholds[threshold_dimension]).astype(int)
	
	return matrix_array

# Binarise per channel
lfp_air_binarised = threshold_binarise_arrays_in_array(lfp_air, channel_medians)

"""
# # Show data for single trial
f, (p1, p2) = plt.subplots(2, sharex=True, sharey=True)
p1.imshow(lfp_air[1], interpolation="none", aspect="auto")
p2.imshow(lfp_air_binarised[1], interpolation="none", cmap="gray", aspect="auto")
plt.show(block=False)
"""

"""
TPM:

Idea is to get phi across time (i.e. to calculate phi at each time bin; how to determine appropriate time bin size? Perception speed?)
Thus, a TPM must be calculated at each time bin (in order to calculate phi)

PHI is calculated using the TPM and state
Method 1:
	To find PHI when awake, find the TPM when awake, use with awake states
	To find PHI when anest, find the TPM when anest, use with anest states
Method 2:
	Find the TPM across conditions (or use awake as baseline, or use anest as baseline?)
	To find PHI when awake, use with awake states
	TO find PHI when anest, use with anest states

Is there some way to quantify difference between TPMs? Maybe awake TPM is very different to anest TPM?
"""

def tpm_window(matrix_array, window_start_sample, window_end_sample):
	"""
	For a given time range (from window_start to window_end, inclusive), build a transitional probability matrix (TPM)
	Iterates through each matrix in the array
	Extracts values at each dimension-1 element, and at each dimension-2 element (corresponding to all channels, and the limited time window, respectively)
	Inputs:
		window_start_sample - which sample (column) to start building from
			If negative, all samples are used (and window_end is ignored)
		window_end_sample - which sample (column) to end building at (it is assumed that window_end doesn't transition to a future state)
	Outputs:
	"""
	
	# Determine number of states
	n_states = 2**matrix_array[0].shape[0] # This assumes binarisation - a node is either ON or OFF
	print(matrix_array[0].shape[0])
	# Declare TPM (all zeros)
	tpm = np.zeros((n_states, n_states))
	
	"""
	TPM Indexing (LOLI):
	e.g. for 4x4 TPM:
	
	0 = 00
	1 = 10
	2 = 01
	3 = 11
	
	Use pyphi.convert.state2loli_index(tuple) to get the index
	
	TODO: Return state-by-node, not state-by-state
	"""
	
	# Declare transition counter (0)
	transition_counter = np.zeros((n_states, 1))
	
	for matrix_counter in range(matrix_array.shape[0]): # For each matrix
		matrix = matrix_array[matrix_counter]
		last_sample = matrix.shape[1] - 1 # Last sample which transitions (i.e. the second last sample)
		
		if window_start_sample == -1:
			window_start = 0
			window_end = last_sample - 1
		else:
			window_start = window_start_sample
			window_end = window_end_sample
		
		for sample_counter in range(window_start, window_end): # For each sample in the matrix (with a future sample)
			sample_current = matrix[:, sample_counter] # Current state (given by the column, assuming that each row is the time series of a channel)
			sample_future = matrix[:, sample_counter+1] # Future state
			
			# Identify current state
			state_current = pyphi.convert.state2loli_index(tuple(np.ndarray.tolist(sample_current)))
			
			# Identify the following state
			state_future = pyphi.convert.state2loli_index(tuple(np.ndarray.tolist(sample_future)))
			
			# Increment TPM transition by 1
			tpm[state_current, state_future] += 1
			
			# Increment transition counter
			transition_counter[state_current] += 1
	
	# Divide elements in TPM by transition counter
	tpm /= transition_counter
	
	return tpm #pyphi.convert.state_by_state2state_by_node(tpm)

"""
# Only use n channels
for trial in range(lfp_air_binarised.size):
	lfp_air_binarised[trial] = lfp_air_binarised[trial][(0, 3), :]

# Build TPM
tpm = tpm_window(lfp_air_binarised, -1, 0)

# Calculate PHI for one trial
trial = lfp_air_binarised[1]
network = pyphi.Network(tpm)
phis = np.zeros((1, trial.shape[1]))
for sample_counter in range(trial.shape[1]):
	print(sample_counter)
	sample = trial[:, sample_counter]
	state = tuple(np.ndarray.tolist(sample))
	subsystem = pyphi.Subsystem(network, state, range(network.size))
	#phis[sample_counter] = pyphi.compute.big_phi(subsystem)
"""

# Only use n channels
for trial in range(lfp_air_binarised.size):
	lfp_air_binarised[trial] = lfp_air_binarised[trial][0:3, :]

# Build TPM
tpm = tpm_window(lfp_air_binarised, -1, 0)

"""
# # Show data for single trial
f, (p1, p2, p3) = plt.subplots(3, sharex=False, sharey=False)
plot1 = p1.imshow(lfp_air[1][0:3, 0:20], interpolation="none", aspect="auto")
#plot1_cbar = f.colorbar(plot1)
plot2 = p2.imshow(lfp_air_binarised[1][:, 0:20], interpolation="none", cmap="gray", aspect="auto")
#plot2_cbar = f.colorbar(plot2)
p3.imshow(tpm, interpolation="none", aspect="auto")
plt.show(block=False)
"""

def independent_joint_tpm(lfp_air_binarised, 0):


f, (p1) = plt.subplots(1, facecolor="white")
plot1 = p1.imshow(tpm, interpolation="none", aspect="auto")
p1.set_xlabel("state")
p1.set_ylabel("state")
plot1_cbar = f.colorbar(plot1)
plot1_cbar.set_label("p")
plt.rcParams.update({'font.size': 22})
plt.show(block=False)

# Calculate PHI for one trial
trial = lfp_air_binarised[1]
network = pyphi.Network(tpm)
phis = np.zeros((1, trial.shape[1]))
for sample_counter in range(trial.shape[1]):
	print(sample_counter)
	sample = trial[:, sample_counter]
	state = tuple(np.ndarray.tolist(sample))
	subsystem = pyphi.Subsystem(network, state, range(network.size))
	#phis[sample_counter] = pyphi.compute.big_phi(subsystem)

def nchannel_phi(matrix_array, n_channels):
	"""
	For each combination of n_channels (n_channels <= number of rows/channels in each matrix), calculates
	the TPM across all trials using all samples, then calculates phi at each transitioning sample

	Assumes that all channels are connected
	
	Inputs:
		matrix_array - array of 2D matrices, all with the same number of elements in the first dimension
		n_channels - how many channels to calculate phi over
	Outputs:
	"""
	
	combo_phis = dict()
	
	# Get the combinations of channels
	channels = np.arange(matrix_array[0].shape[0]) # Get list of channels
	channel_combos = [list(combination) for combination in itertools.combinations(channels, n_channels)] # Get combos of channels
	
	for combo in channel_combos:
		print(combo)
		# Extract relevant channels (i.e. channels specificed by the combo) across all trials
		samples_all = matrix_array[0][combo, :]
		for matrix_counter in range(1, matrix_array.size):
			np.concatenate((samples_all, matrix_array[matrix_counter][combo, :]), axis=1)
		# Build the TPM from those channels
		tpm_input = np.empty((1), dtype=object)
		tpm_input[0] = samples_all
		tpm = tpm_window(tpm_input, -1, 0) # We place the concatenated matrix into an array as input for tpm_window(), which takes an array of matrices
		
		print("computing phi")
		# Compute phi at each trial sample
		network = pyphi.Network(tpm)
		phis = np.empty((matrix_array.size), dtype=object)
		for trial_counter in range(matrix_array.size):
			phis[trial_counter] = np.empty(matrix_array[trial_counter].shape[1])
			trial = matrix_array[trial_counter]
			for sample_counter in range(matrix_array[trial].shape[1]):
				print("trial" + str(trial_counter) + " sample" + str(sample_counter))
				sample = trial[combo, sample_counter]
				print(sample)
				state = tuple(np.ndarray.tolist(sample))
				print("state done")
				subsystem = pyphi.Subsystem(network, state, range(network.size))
				print("subsystem done")
				phis[trial_counter][sample_counter] = pyphi.compute.big_phi(subsystem)
				print("phi = " + str(phis[trial_counter][sample_counter]))
		
		# Add combo results to dictionary
		combo_phis[combo] = phis
	
	return combo_phis

"""
Conditions at which to compute phi:
View all time and trials consecutively as a line graph - one matrix
View all time and trials as an image/matrix, trim as required before viewing - one matrix per trial
	If using this, then you can reconstruct the first option
Average across trials (time is the length of a trial, trim as required before averaging) - one matrix per trial

"""

def recursive_phi(max_channels):
	"""
	Computes phi over the power set of channels, from using 2 channels to up to max_channels
	"""
	return 1

# phis = nchannel_phi(lfp_air_binarised, 2)

def plot_mean(matrix):
	"""
	Calculates the mean and standard deviation across the first dimension (axis 0), and plots across the second dimension (axis 1)
	Inputs: numpy 2D array
	"""
	if len(matrix.shape) > 2:
		return
	
	if len(matrix.shape) == 2: # two dimensions provided
		x_values = np.arange(0, matrix.shape[1])
		mean_values = matrix.mean(axis=0)
		std_values = matrix.std(axis=0)
	else: # only one dimension was provided
		x_values = np.arange(0, matrix.shape[0])
		mean_values = matrix
		std_values = matrix
	
	plt.plot(x_values, mean_values)
	if len(matrix.shape) == 2: # two dimensions provided
		plt.fill_between(x_values, mean_values - std_values, mean_values + std_values, alpha=0.5, linewidth=0)
	plt.autoscale(enable=True, axis='x', tight=True)
	plt.show(block=False)

def main():
	print("main() running")
	lfp = fly_data["LFP"][0, 0][0, :]
	#plot_mean(lfp)
	print("main() finished")


if __name__ == '__main__':
	main()

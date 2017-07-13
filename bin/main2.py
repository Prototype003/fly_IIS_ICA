"""
To find package install locations, use:
import site; site.getsitepackages()

Anaconda (not standardalone Python):
Location should be "C:/Users/this_/Anaconda3/Lib/site-packages"
pyphi is installed in the above directory
"""

print("Setting up environment")

from fly_phi import load_mat, save_mat, concat_arrays_in_array, threshold_binarise_arrays_in_array, tpm_window, nchannel_phi, plot_mean
import numpy as np
import pyphi

print("packages imported")

def main():
	data_directory = "../../flies/fly_data/trials_anesthDescAdded11092014_bPlrRerefTyp2/Analyzed_WalkingDrosDror115610072014/" # 8 channels (after bipolar re-referencing - need confirmation, should create function for re-referencing)
	#data_directory = "../../flies/fly_data/trials_anesthDescAdded11092014/Analyzed_WalkingDrosDror115610072014/" # 16 channels (raw)
	
	loaded_fly_data = load_mat(data_directory)
	fly_data = loaded_fly_data["trials"]
	print("Data format: " + str(type(fly_data)) + " : " + str(fly_data.shape))

	lfp = fly_data["LFP"][0, :] # LFP is a 1D vector (but numpy has it as a 2D matrix with the size of the first dimension=1) [0, :] gets the first (in this case, only) element in the 1st dimension, and all in the 2nd dimension
	anesth_label = fly_data["anesthLabel"][0, :]
	on_trial = fly_data["onTrial"][0, :]


	# Get all rests (times between on trials) for the air condition
	lfp_air = lfp[np.ix_(np.logical_and(anesth_label==1, on_trial==0))]


	# Concatenate all trials together (for finding the global median - for binarisation)
	lfp_air_concatenated = concat_arrays_in_array(lfp_air)


	# Find global median per channel (using concatenated trials)
	channel_medians = np.median(lfp_air_concatenated, axis=1) # Median across samples, per channel


	# Binarise per channel
	lfp_air_binarised = threshold_binarise_arrays_in_array(lfp_air, channel_medians)


	# Only use n channels
	for trial in range(lfp_air_binarised.size):
		lfp_air_binarised[trial] = lfp_air_binarised[trial][(0, 3), :]

	# Build TPM
	tpm = tpm_window(lfp_air_binarised, -1, 0)

	# Calculate PHI for one trial
	# trial = lfp_air_binarised[1]
	# network = pyphi.Network(tpm)
	# phis = np.zeros((trial.shape[1]))
	# for sample_counter in range(trial.shape[1]):
		# print(sample_counter)
		# sample = trial[:, sample_counter]
		# state = tuple(np.ndarray.tolist(sample))
		# subsystem = pyphi.Subsystem(network, state, range(network.size))
		# phis[sample_counter] = pyphi.compute.big_phi(subsystem)

	phis = nchannel_phi(lfp_air_binarised, 2)
	
	results_directory = "../results/phi_results.mat"
	results = dict()
	results["phis"] = phis
	save_mat(results_directory, results)


if __name__ == '__main__':
	main()

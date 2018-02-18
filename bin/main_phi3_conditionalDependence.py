'''
This script is for demonstrating the issue of conditional independence regarding Phi3 and PyPhi
'''

setSize_min = 2


# Imports
import sys
import pyphi as pyphi
import numpy as np
import scipy.io as sio
import math as math
from fly_phi import *

# Load phi results
data_directory = "results/"
data_file = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree.mat"

results_directory = "analysis_results/"
results_file = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree_tpmDependence.mat"

loaded_data = load_mat(data_directory + data_file)
phis = loaded_data['phis']
print("Phis loaded")
sys.stdout.flush()

# This gives us an indication of which TPMs are changed more after conversion and back-conversion
def tpm_comparison(tpm1, tpm2):
	# Assumes that the sizes of tpm1 and tpm2 are the same
	# Returns:
	#	Sum of absolute differences between all probabilities
	#	Pearson correlation coefficient between the two tpms
	
	# Flatten tpms into 1D
	tpm1 = tpm1.flatten()
	tpm2 = tpm2.flatten()
	
	diff = 0
	for prob in range(0, tpm1.size):
		diff = diff + abs(tpm1[prob] - tpm2[prob])
	
	corr = np.corrcoef(tpm1, tpm2)
	corr = corr[0, 1]
	
	return [diff, corr]

# Results storage structure
tpm_results = [dict() for n in range(phis.size)]
for setSize_counter in range(0, phis.size):
	nChannels = setSize_counter + setSize_min
	tpm_results[setSize_counter]['sbs1'] = np.zeros(phis[0, setSize_counter][0, 0]['tpms'].shape)
	tpm_results[setSize_counter]['sbn'] = np.zeros((2**nChannels, nChannels) + phis[0, setSize_counter][0, 0]['tpms'].shape[2:])
	tpm_results[setSize_counter]['sbs2'] = np.zeros(phis[0, setSize_counter][0, 0]['tpms'].shape)
	tpm_results[setSize_counter]['diff'] = np.zeros(phis[0, setSize_counter][0, 0]['tpms'].shape[2:])
	tpm_results[setSize_counter]['corr'] = np.zeros(phis[0, setSize_counter][0, 0]['tpms'].shape[2:])

# Iterate through all TPMs and find how much they change after conversion and back-conversion
for setSize_counter in range(0, phis.size):
	nChannels = setSize_counter + setSize_min
	for set_counter in range(0, phis[0, setSize_counter][0, 0]['tpms'].shape[2]):
		for fly_counter in range(0, phis[0, setSize_counter][0, 0]['tpms'].shape[3]):
			for condition_counter in range(0, phis[0, setSize_counter][0, 0]['tpms'].shape[4]):
				for tau_counter in range(0, phis[0, setSize_counter][0, 0]['tpms'].shape[5]):
					print("nChannels" + str(nChannels) + "set" + str(set_counter) + "fly" + str(fly_counter) + "condition" + str(condition_counter) + "tau" + str(tau_counter))
					sys.stdout.flush()
					
					# Conversions
					sbs1 = phis[0, setSize_counter][0, 0]['tpms'][:, :, set_counter, fly_counter, condition_counter, tau_counter]
					sbn = pyphi.convert.state_by_state2state_by_node(sbs1)
					sbs2 = pyphi.convert.state_by_node2state_by_state(sbn)
					
					# Find difference between original and back-converted tpms
					comparison = tpm_comparison(sbs1, sbs2)
					
					# Store results
					tpm_results[setSize_counter]['sbs1'][:, :, set_counter, fly_counter, condition_counter, tau_counter] = sbs1
					tpm_results[setSize_counter]['sbn'][:, :, set_counter, fly_counter, condition_counter, tau_counter] = sbn.reshape((2**nChannels, nChannels), order='F')
					tpm_results[setSize_counter]['sbs2'][:, :, set_counter, fly_counter, condition_counter, tau_counter] = sbs2
					tpm_results[setSize_counter]['diff'][set_counter, fly_counter, condition_counter, tau_counter] = comparison[0]
					tpm_results[setSize_counter]['corr'][set_counter, fly_counter, condition_counter, tau_counter] = comparison[1]

# Save converted and backconverted TPMs
save_mat(results_directory + results_file, {'tpms': tpm_results})
print('saved')
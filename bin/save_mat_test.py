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
flies = np.arange(0, 13)
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
	"_phithree"

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
	
	# Initialise phi_matrix (channel sets x trials x flies x conditions x taus x nBins)
	phi_threes = np.zeros((len(channel_sets), fly_data.shape[2], fly_data.shape[3], fly_data.shape[4], len(taus), len(nBins)))
	
	# Storage of state counts (states x channel sets x trials x flies x conditions x taus x nBins) and phis (-trials and -nBins dimensions)
	state_counter = np.zeros((n_states, len(channel_sets), fly_data.shape[2], fly_data.shape[3], fly_data.shape[4], len(taus), len(nBins)))
	state_phis = np.zeros((n_states, len(channel_sets), fly_data.shape[3], fly_data.shape[4], len(taus), len(nBins))) # This will hold the phi values of each state
	
	# Storage of state trajectory (samples x channel sets x trials x flies x conditions x taus x nBins)
	#state_trajectories = [np.zeros((fly_data.shape[0], fly_data.shape[2], fly_data.shape[3], fly_data.shape[4], len(taus), len(nBins))) for n in range(len(channel_sets))]
	
	# Storage of TPMs (states x states x channel sets x flies x conditions x taus)
	tpms = np.zeros((n_states, n_states, len(channel_sets), fly_data.shape[3], fly_data.shape[4], len(taus)))

	phis[nChannels_counter]['phi_threes'] = phi_threes
	phis[nChannels_counter]['state_counters'] = state_counter
	phis[nChannels_counter]['state_phis'] = state_phis
	#phis[nChannels_counter]['state_trajectories'] = state_trajectories
	phis[nChannels_counter]['tpms'] = tpms

# Save ###########################################################################

print('saving')
"""
for nChannels_counter in range(0, len(nChannels)):
	results_file = data_file + \
	"_detrend" + str(prep_detrend) + \
	"_zscore" + str(prep_zscore) + \
	"_nChannels" + str(nChannels[nChannels_counter]) + \
	"_phithree"
	save_mat(results_directory+results_file, \
		{'nChannels': phis[nChannels_counter]['nChannels'], \
		'channel_sets': phis[nChannels_counter]['channel_sets'], \
		'taus': phis[nChannels_counter]['taus'], \
		'nBins': phis[nChannels_counter]['nBins'], \
		'phi_threes': phis[nChannels_counter]['phi_threes'], \
		'state_counters': phis[nChannels_counter]['state_counters'], \
		'state_phis': phis[nChannels_counter]['state_phis'], \
		#'state_trajectories': phis[nChannels_counter]['state_trajectories'], \
		'tpms': phis[nChannels_counter]['tpms'] \
		})
"""
save_mat(results_directory+results_file, {'phis': phis})

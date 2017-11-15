
# This script is to concatenate results across flies

import sys
import numpy as np
from fly_phi import *

prefix_common = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_allPartitions'
infix_common = 'fly'
suffix_common = '.mat'

data_directory = 'results/preformatted_results/'

results_directory = 'results/'
results_filename = prefix_common + suffix_common

flies = np.arange(0, 13)

# Load first fly - all other fly results will be concatenated to this fly
print(flies[0])
sys.stdout.flush()
phis = load_mat(data_directory + prefix_common + infix_common + str(flies[0]) + suffix_common)
phis = phis['phis']
phis_data = phis[0]

#sys.exit()

# Concatenate other fly results
for fly in range(1, len(flies)): # Start from the second fly
	print(fly)
	sys.stdout.flush()
	
	# Load fly data
	fly_data = load_mat(data_directory + prefix_common + infix_common + str(flies[fly]) + suffix_common)
	fly_data = fly_data['phis']
	fly_data = fly_data[0]
	
	# Concatenate results to previous flies
	phis_data[0][0, 0]['phi_threes'] = np.concatenate((phis_data[0][0, 0]['phi_threes'], fly_data[0][0, 0]['phi_threes']), axis=2)
	phis_data[0][0, 0]['mips'] = np.concatenate((phis_data[0][0, 0]['mips'], fly_data[0][0, 0]['mips']), axis=2)
	phis_data[0][0, 0]['state_counters'] = np.concatenate((phis_data[0][0, 0]['state_counters'], fly_data[0][0, 0]['state_counters']), axis=3)
	phis_data[0][0, 0]['state_phis'] = np.concatenate((phis_data[0][0, 0]['state_phis'], fly_data[0][0, 0]['state_phis']), axis=2)
	phis_data[0][0, 0]['tpms'] = np.concatenate((phis_data[0][0, 0]['tpms'], fly_data[0][0, 0]['tpms']), axis=3)
	phis_data[0][0, 0]['state_partitions'] = np.concatenate((phis_data[0][0, 0]['state_partitions'], fly_data[0][0, 0]['state_partitions']), axis=3)
	phis_data[0][0, 0]['state_partitions_phis'] = np.concatenate((phis_data[0][0, 0]['state_partitions_phis'], fly_data[0][0, 0]['state_partitions_phis']), axis=3)
	
phis[0] = phis_data

# Save concatenated results
print('saving')
sys.stdout.flush()
save_mat(results_directory+results_filename, {'phis': phis})
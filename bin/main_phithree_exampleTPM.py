# This is for creating a figure of a worked example for a simple system of 2 channels
# Take a channel pair from a fly which gives high phi
# Present the empirical TPM, and rebuild an independent TPM from the actual channel data (by multiplying independent channel probabilities)
# Present the associated phi values from using each TPM

import matplotlib.pyplot as plt
import sklearn.preprocessing as skp
from fly_phi import *



# Parameter Selection ################################################
# Channel pair which gives max delta phi in fly 1:
#	Channel set 51 gives:
#		delta = 0.0076
#		air = 0.0102
#		iso = 0.0025
#	Channels (0-indexed) = [4 5]

fly = 1
channel_set = [4, 5]
trial = 1 # Don't know which trial is best
condition = 1 # air condition should have integration
tau = 4

samples = np.arange(0, 20)
######################################################################

# Get Data ###########################################################
data_directory = "workspace_results/"
data_file = "split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim"

print("Loading fly data")
loaded_data = load_mat(data_directory + data_file + ".mat")
fly_data = loaded_data['fly_data'] # Reminder: dimensions are: (samples x channels x trials x flies x conditions)
print("Fly data loaded")
######################################################################

# Filter for data from 1 fly only
#channel_data = fly_data[:, channel_set, 1, fly, :]

# Discretise
fly_data_discretised, n_values, channel_medians, *unused = binarise_global_median(fly_data)
print("Fly data discretised")

# Filter for data from 1 fly only - preferably in the air condition
channel_data = fly_data_discretised[:, :, trial, fly, condition]
channel_data = channel_data[:, channel_set]
channel_data = channel_data[samples, :] # I split these lines from [:, channel_set, trial, fly, condition] because specifying all the indices in one line reorders the dimensions in the output for whatever reason

# Build TPM from whole system
tpm = build_tpm(channel_data, tau, n_values)

# Build TPM per channel

# Multiple per-channel TPMs to obtain TPM for cut system


# Plot data
f, (p1) = plt.subplots(1, facecolor="white")
plot1 = p1.imshow(np.transpose(channel_data[:, :, 1]), cmap='gray', interpolation="none", aspect="auto")
plt.show(block=False)
# Plot real TPM
f, (p1) = plt.subplots(1, facecolor="white")
plot1 = p1.imshow(tpm, interpolation="none", aspect="auto")
p1.set_xlabel("state")
p1.set_ylabel("state")
plot1_cbar = f.colorbar(plot1)
plot1_cbar.set_label("p")
plt.rcParams.update({'font.size': 22})
plt.show(block=False)
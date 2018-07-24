function [mip, phi] = compute_discrete_nonGlobal(phi_type, fly, condition, trial, nChannels, set_counter, tau, save_file)
% Assumes data is already discretised to states, starting from 0 (e.g. 0, 1)
%
% compute_discrete_nonGlobal(phi_type, fly, condition, trial, nChannels, set_counter, tau, save_file)
%
% Inputs:
%   phi_type:
%       'SI': phi_H, stochastic interaction
%       'star': phi_star, based on mismatched decoding
%       Use 'help MIP_search' for all options

file_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';
file_suffix = '.mat';

% Add paths
bin_location = '../';
addpath(genpath([bin_location 'PhiToolbox-master/']));
data_dir = [bin_location 'fly_data_split/'];
results_dir = ['results_split/'];
results_file = [file_prefix...
    '_nChannels' num2str(nChannels)...
    '_phi' phi_type...
    '_f' num2str(fly, '%02.f') 'c' num2str(condition) 'tau' num2str(tau) 's' num2str(set_counter, '%04.f') 't' num2str(trial)...
    file_suffix];

% Load data file
data_file = [file_prefix...
    '_f' num2str(fly) 'c' num2str(condition) 't' num2str(trial)...
    file_suffix];
loaded = load([data_dir data_file]);

% Determine channels to use
channel_sets = nchoosek((1:size(loaded.data, 2)), nChannels);
channel_set = channel_sets(set_counter, :);

% Get data from channels
data = loaded.data(:, channel_set)' + 1; % channels x samples; increment all states to start from 1

% Phi-star options
options.type_of_dist = 'discrete';
options.type_of_phi = phi_type;
options.type_of_MIPsearch = 'Exhaustive';

% Parameters
params.tau = tau;
params.number_of_states = max(data(:));

[mip_old, phi] = MIP_search(data, params, options);
% mip_old format: each index corresponds to a channel, the element gives the partition which the channel belongs to, in the MIP

% Reformat mip
mip = partition_z2cell(mip_old);

% Save
if save_file == 1
    save([results_dir results_file], 'phi', 'mip');
end

end


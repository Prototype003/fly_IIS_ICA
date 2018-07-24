function [mip, phi] = compute_discrete_global(phi_type, fly, condition, nChannels, set_counter, tau, save_file)
% Takes the raw data and binarises it before computing phi
%
% compute_discrete_global(phi_type, fly, condition, nChannels, set_counter, tau, save_file)
%
% Inputs:
%   phi_type:
%       'SI': phi_H, stochastic interaction
%       'star': phi_star, based on mismatched decoding
%       Use 'help MIP_search' for all options
%   tau:
%       Integer - actual tau value
%   save_file:
%       Integer - if 1, saves results to file in directory 'results_split/'

file_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';
file_suffix = '.mat';

% Add paths
bin_location = '../';
addpath(genpath([bin_location 'PhiToolbox-master/']));
data_dir = [bin_location 'workspace_results/'];
results_dir = ['results_split/'];
results_file = [file_prefix...
    '_nChannels' num2str(nChannels)...
    '_phi' phi_type...
    '_f' num2str(fly, '%02.f') 'c' num2str(condition) 'tau' num2str(tau) 's' num2str(set_counter, '%04.f')...
    file_suffix];

% Load data file
data_file = [file_prefix file_suffix];
loaded = load([data_dir data_file]);

% Determine channels to use
channel_sets = nchoosek((1:size(loaded.fly_data, 2)), nChannels);
channel_set = channel_sets(set_counter, :);

% Get data from channels
data = loaded.fly_data(:, channel_set, :, fly, condition); % samples x channels x trials
data = permute(data, [2 1 3]); % channels x samples x trials

% Join trials
data_joined = zeros(size(data, 1), size(data, 2) * size(data, 3)); % channels x samples
sample_counter = (1 : size(data, 2));
for trial = 1 : size(data, 3)
    data_joined(:, sample_counter) = data(:, :, trial);
    sample_counter = sample_counter + size(data, 2);
end

% Binarise data based on median split
medians = median(data_joined, 2);
medians = repmat(medians, [1 size(data_joined, 2)]);
data_joined = data_joined > medians;
data = data_joined + 1; % increment states so that lowest state is 1

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


%% DESCRIPTION

%{

Get channels whose average 2ch phi predict the highest 4ch phi, for 1 fly

%}

%% Setup

measure = 'phi_three';
tau = 1; % 1 = 4ms; 2 = 8ms; 3 = 16ms
if tau == 1
    tau_string = '4';
elseif tau == 2
    tau_string = '8';
elseif tau == 3
    tau_string = '16';
end

freq_range = (1:42); %(1:83); % corresponding to ~5Hz and ~10Hz, check the 'frequencies' vector
freq_range_string = '0-5Hz'; %'0-10Hz';

fontsize = 11; % Used for drawing label letters

bin_location = '../';
addpath(bin_location);
addpath([bin_location 'figure_code/']);

results_directory = [bin_location 'workspace_results/'];

%% Load phi values (WITHIN)

% We'll not rename this because the filesize is big and previous figures use this as well
[phis, measure_string] = phi_load(measure, 0, bin_location);

%% Load phi values (ACROSS)

[phis_a, measure_string] = phi_load(measure, 1, bin_location);

%% Get mean values per channel

channels = (1:15);
channel_means = cell(length(phis), 1); % We will save channel means so we can use them for selecting large sets of channels

% Storage of channel means
for nChannels_counter = 1 : length(phis)
    channel_sets = phis{nChannels_counter}.channel_sets;
    
    % Phi values
    tmp_size = size(phis{nChannels_counter}.phis(:, :, :, :, tau));
    tmp_size(1) = length(channels);
    phis{nChannels_counter}.channel_sums = zeros(tmp_size);
    phis{nChannels_counter}.channel_means = zeros(tmp_size);
    phis{nChannels_counter}.set_counters = zeros(tmp_size);
    phis_a{nChannels_counter}.channel_sums = zeros(tmp_size);
    phis_a{nChannels_counter}.channel_means = zeros(tmp_size);
    phis_a{nChannels_counter}.set_counters = zeros(tmp_size);
    
    % Sum values for each channel (sum across networks which contain the channel)
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                
                phis{nChannels_counter}.channel_sums(channel, :, :, :) = phis{nChannels_counter}.channel_sums(channel, :, :, :) + phis{nChannels_counter}.phis(channel_set, :, :, :, tau);
                phis{nChannels_counter}.set_counters(channel, :, :, :) = phis{nChannels_counter}.set_counters(channel, :, :, :) + 1;
                
                phis_a{nChannels_counter}.channel_sums(channel, :, :, :) = phis_a{nChannels_counter}.channel_sums(channel, :, :, :) + phis_a{nChannels_counter}.phis(channel_set, :, :, :, tau);
                phis_a{nChannels_counter}.set_counters(channel, :, :, :) = phis_a{nChannels_counter}.set_counters(channel, :, :, :) + 1;
                
            end
        end
    end
    
    % Average values
    phis{nChannels_counter}.channel_means = phis{nChannels_counter}.channel_sums ./ phis{nChannels_counter}.set_counters;
    phis_a{nChannels_counter}.channel_means = phis_a{nChannels_counter}.channel_sums ./ phis_a{nChannels_counter}.set_counters;
    
end

%%

fly = 1;
condition = 1;
tau = 1;
nChannels_counter = 1;

values = phis_a{nChannels_counter}.phis(:, :, fly, condition, tau);
channel_means = mean(phis_a{1}.channel_means(:, :, fly, condition), 2); % values are identical across trials (because they're a duplicate of the same one trial)
channel_sets = phis_a{nChannels_counter}.channel_sets;

% Build values from average of single channel average values
values_rebuilt = zeros(size(values));
for channel_set = 1 : size(channel_sets, 1)
    channels = channel_sets(channel_set, :);
    values_rebuilt(channel_set) = mean(channel_means(channels));
end

[phi_max, index] = max(values_rebuilt);

channel_set = channel_sets(index, :)

figure;
scatter(values_rebuilt, values, '.');

%% Calculate phi-star

addpath('../'); % We are using the phistar_mip.m function

% Load data
data_directory = '../workspace_results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';
load([data_directory data_file '.mat']);

tau = 4;
% condition and fly are the same as in previous section

% Flatten trials
data_flat = fly_data(:, :, 1, :, :);
for trial = 2 : size(fly_data, 3)
    data_flat = cat(1, data_flat, fly_data(:, :, trial, :, :));
end
fly_data = data_flat;

% Compute covariances
[cov_past_past, cov_past_present, cov_present_present] = Cov_comp(fly_data(:, channel_set, :, fly, condition)', tau);

[phi_tmp,...
    h,...
    hcond,...
    mi,...
    mi_star,...
    mip,...
    partitions,...
    partitions_phi_tmp,...
    partitions_hs,...
    partitions_hconds,...
    partitions_mis,...
    partitions_mi_stars...
    ] = phistar_mip(cov_past_past, cov_past_present, cov_present_present, channel_set);

% phi_stars(channel_set, trial, fly, condition, tau_counter, nBins_counter) = phi_tmp(1);
% phi_stars_normalised(channel_set, trial, fly, condition, tau_counter, nBins_counter) = phi_tmp(2);
% 
% partitions_phi_stars(:, channel_set, trial, fly, condition, tau_counter, nBins_counter) = partitions_phi_tmp(:, 1);
% partitions_phi_stars_normalised(:, channel_set, trial, fly, condition, tau_counter, nBins_counter) = partitions_phi_tmp(:, 2);


%% DESCRIPTION

%{
This script computes phi_star for a combinations of n channels, for some
set of taus.

Note: it only supports 1 bin per epoch, i.e. the bin length = epoch length
- tau
%}

%% SETUP

prep_detrend = 0;
prep_zscore = 0;
prep_medianSplit = 0;

flies = (1:13); % This will be parallelised
nChannels = (2:4); % This determines how many channels to consider at a time
nBins = [1]; % Script currently only supports the case of 1 bin (i.e. epoch length - tau)
taus = [4, 8, 16];

nFlies = length(flies);

data_directory = 'workspace_results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

results_directory = 'results/';
results_filename = [data_file '_detrend' num2str(prep_detrend) '_zscore' num2str(prep_zscore) '_nChannels' num2str(nChannels(1)) 't' num2str(nChannels(end)) '_medianSplit' num2str(prep_medianSplit) '_phistar_allPartitions']; %flies' num2str(flies(1)) 't' num2str(flies(end))];

%% LOAD

load([data_directory '/' data_file '.mat']);

%% START

if prep_detrend == 1
    for condition = 1 : 2
        for fly = 1 : size(fly_data, 4)
            for trial = 1 : size(fly_data, 3)
                fly_data(:, :, trial, fly, condition) = detrend(fly_data(:, :, trial, fly, condition));
            end
        end
    end
end
if prep_zscore == 1
    fly_data = zscore(fly_data, [], 1);
end
if prep_medianSplit == 1
    fly_data = binarise_global_median(fly_data);
end

%{
Storage structure:

Cell array, each cell holds results for nChannels (e.g. 2 channels, 3,
channels, etc)

Struct should contain:
    channel_sets - holds all channel sets (matrix: combo x nChannels)
    taus - holds the tau lags (vector: taus)
    nBins - holds the number of bins (vector: nBins)
    nChannels - number of channels used in computation

    phis - holds the phi values (matrix: channel sets x trials x flies x conditions x taus x nBins)
%}
phis = cell(length(nChannels), 1);

% Compute phi per channel combination


for nChannels_counter = 1 : length(nChannels)
    channels = nChannels(nChannels_counter);
    channel_sets = nchoosek((1:size(fly_data, 2)), channels);
    
    % Storage
    % channels sets x trials x flies x conditions x taus x nBins
    phi_stars = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    phi_stars_normalised = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    mis = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    hconds = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    hs = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    mi_stars = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    mips = cell(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    
    % Storage for every partition
    partitions = cell(Bell(channels)-1, size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    partitions_phi_stars = zeros(Bell(channels)-1, size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    partitions_phi_stars_normalised = zeros(Bell(channels)-1, size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    partitions_mis = zeros(Bell(channels)-1, size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    partitions_mi_stars = zeros(Bell(channels)-1, size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    partitions_hs = zeros(Bell(channels)-1, size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    partitions_hconds = zeros(Bell(channels)-1, size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    
    for fly_counter = 1 : nFlies
        fly = flies(fly_counter);
        disp(['nChannels' num2str(channels) ' fly' num2str(fly)]); 
        for nBins_counter = 1 : length(nBins)
            for tau_counter = 1 : length(taus)
                tau = taus(tau_counter);
                for channel_set = 1 : size(channel_sets, 1)
                    disp(['Channel Set ' num2str(channel_set)]);
                    for condition = 1 : size(fly_data, 5)
                        for trial = 1 : size(fly_data, 3)
                            %disp(['fly' num2str(fly) ' tau' num2str(taus(tau_counter)), ' combo' num2str(channel_set) ' condition' num2str(condition) ' trial' num2str(trial)]); % Keep in mind this is slow
                            % Compute covariances
                            [cov_present_present, cov_present_past, cov_past_past] = Cov_comp_sample(fly_data(:, channel_sets(channel_set, :), trial, fly, condition)', tau);

                            % Compute phi (find the MIP)
                            [phi_tmp,...
                                hs(channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                hconds(channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                mis(channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                mi_stars(channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                mips{channel_set, trial, fly, condition, tau_counter, nBins_counter},...
                                partitions(:, channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                partitions_phi_tmp,...
                                partitions_mis(:, channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                partitions_mi_stars(:, channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                partitions_hs(:, channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                partitions_hconds(:, channel_set, trial, fly, condition, tau_counter, nBins_counter)...
                                ] = phistar_mip(cov_past_past, cov_present_past, cov_present_present, channel_sets(channel_set, :));
                            
                            phi_stars(channel_set, trial, fly, condition, tau_counter, nBins_counter) = phi_tmp(1);
                            phi_stars_normalised(channel_set, trial, fly, condition, tau_counter, nBins_counter) = phi_tmp(2);
                            
                            partitions_phi_stars(:, channel_set, trial, fly, condition, tau_counter, nBins_counter) = partitions_phi_tmp(:, 1);
                            partitions_phi_stars_normalised(:, channel_set, trial, fly, condition, tau_counter, nBins_counter) = partitions_phi_tmp(:, 2);
                            
                        end
                    end
                end
            end
        end
    end
    
    phis{nChannels_counter} = struct();
    phis{nChannels_counter}.phi_stars = phi_stars;
    phis{nChannels_counter}.phi_stars_normalised = phi_stars_normalised;
    phis{nChannels_counter}.mis = mis;
    phis{nChannels_counter}.hconds = hconds;
    phis{nChannels_counter}.hs = hs;
    phis{nChannels_counter}.mi_stars = mi_stars;
    phis{nChannels_counter}.mips = mips;
    
    phis{nChannels_counter}.partitions = partitions;
    phis{nChannels_counter}.partitions_phi_stars = partitions_phi_stars;
    phis{nChannels_counter}.partitions_phi_stars_normalised = partitions_phi_stars_normalised;
    phis{nChannels_counter}.partitions_mis = partitions_mis;
    phis{nChannels_counter}.partitions_mi_stars = partitions_mi_stars;
    phis{nChannels_counter}.partitions_hs = partitions_hs;
    phis{nChannels_counter}.partitions_hconds = partitions_hconds;
    
    phis{nChannels_counter}.nChannels = channels;
    phis{nChannels_counter}.channel_sets = channel_sets;
    phis{nChannels_counter}.taus = taus;
    phis{nChannels_counter}.nBins = nBins;
    
end

%% SAVE

if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename '.mat'], 'phis', '-v7.3');

% %% Function: global median split
% 
% function [binarised] = binarise_global_median(fly_data)
% % MATLAB conversion of binarise_global_median() in fly_phi.py
% % Finds the median for each channel, across all epoch-trials
% % Binarises samples based on the median: 1 if greater than median, 0 otherwise
% 
% % Reorder dimensions to prepare for reshape
% reordered = permute(fly_data, [1 3 2 4 5]);
% % Reformat - collapse on trials so that the first dimension is samples x trials
% sizes = size(reordered);
% globalised = reshape(reordered, [sizes(1)*sizes(2) sizes(3) sizes(4) sizes(5)]);
% 
% % Get median per channel, fly, condition
% channel_medians = median(globalised, 1);
% 
% % Binarise based on median:
% % Greater than median = 1
% % Less than or equal to median = 0
% binarised = globalised > channel_medians;
% 
% % Return to original shape
% binarised = reshape(binarised, [sizes(1) sizes(2) sizes(3) sizes(4) sizes(5)]);
% binarised = permute(binarised, [1 3 2 4 5]);
% 
% end
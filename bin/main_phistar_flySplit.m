%% Description

%{
Splits large phis structure into several - 1 per fly
%}

%% Setup

data_directory = 'results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_medianSplit0_phistar_allPartitions.mat';

prefix_common = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_medianSplit0_phistar_allPartitions';
infix = 'fly';
suffix_common = '.mat';

results_directory = 'results/';

%% Load

load([data_directory data_file]);

%% Get basic parameters
% Assumes number of flies is constant across all nChannels 

nFlies = size(phis{1}.mips, 3);

%% Split into separate flies

phis_all = phis;

for fly = 1 : nFlies
    disp(fly);
    
    phis = cell(1, length(phis_all));
    
    results_file = [prefix_common infix num2str(fly) suffix_common];
    
    for nChannels_counter = 1 : length(phis_all)
        phis{nChannels_counter} = struct();
        
        % Get constants
        nChannels = phis_all{nChannels_counter}.nChannels;
        channel_sets = phis_all{nChannels_counter}.channel_sets;
        taus = phis_all{nChannels_counter}.taus;
        nBins = phis_all{nChannels_counter}.nBins;

        % Extract all data for this fly
        phis{nChannels_counter}.phi_stars = phis_all{nChannels_counter}.phi_stars(:, :, fly, :, :, :);
        phis{nChannels_counter}.phi_stars_normalised = phis_all{nChannels_counter}.phi_stars_normalised(:, :, fly, :, :, :);
        phis{nChannels_counter}.mis = phis_all{nChannels_counter}.mis(:, :, fly, :, :, :);
        phis{nChannels_counter}.hconds = phis_all{nChannels_counter}.hconds(:, :, fly, :, :, :);
        phis{nChannels_counter}.hs = phis_all{nChannels_counter}.hs(:, :, fly, :, :, :);
        phis{nChannels_counter}.mi_stars = phis_all{nChannels_counter}.mi_stars(:, :, fly, :, :, :);
        phis{nChannels_counter}.mips = phis_all{nChannels_counter}.mips(:, :, fly, :, :, :);
        
        phis{nChannels_counter}.partitions = phis_all{nChannels_counter}.partitions(:, :, :, fly, :, :, :);
        phis{nChannels_counter}.partitions_phi_stars = phis_all{nChannels_counter}.partitions_phi_stars(:, :, :, fly, :, :, :);
        phis{nChannels_counter}.partitions_phi_stars_normalised = phis_all{nChannels_counter}.partitions_phi_stars_normalised(:, :, :, fly, :, :, :);
        phis{nChannels_counter}.partitions_mis = phis_all{nChannels_counter}.partitions_mis(:, :, :, fly, :, :, :);
        phis{nChannels_counter}.partitions_mi_stars = phis_all{nChannels_counter}.partitions_mi_stars(:, :, :, fly, :, :, :);
        phis{nChannels_counter}.partitions_hs = phis_all{nChannels_counter}.partitions_hs(:, :, :, fly, :, :, :);
        phis{nChannels_counter}.partitions_hconds = phis_all{nChannels_counter}.partitions_hconds(:, :, :, fly, :, :, :);
        
        phis{nChannels_counter}.nChannels = phis_all{nChannels_counter}.nChannels;
        phis{nChannels_counter}.channel_sets = phis_all{nChannels_counter}.channel_sets;
        phis{nChannels_counter}.taus = phis_all{nChannels_counter}.taus;
        phis{nChannels_counter}.nBins = phis_all{nChannels_counter}.nBins;
    end
    
    % Save
    save([results_directory results_file], 'phis');
    
end



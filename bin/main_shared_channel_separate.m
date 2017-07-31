%% Description

%{
Separate raw phi results into two sections:
    Channel sets in which no pair of channels shared a bipolar re-referencing derivation
    Channel sets in which at least one pair of channels shared a bipolar re-referencing derivation

e.g. 1 and 3 do not share a common electrode, but 1 and 2 do

%}

%% Setup

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    ];

% CAREFUL! Don't want to overwrite the original results, so ensure this name is different
% data_filename!!! (_shareFiltered)
results_directory = 'results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_shareFiltered'
    ];

%% Load

disp('loading');

% Phi-3
% .mips (holds the cuts) dimensions: states x sets x flies x conditions x taus
load([data_directory data_filename '_phithree.mat']);
phi_threes = phis;

% Phi-star
% .mips (holds the mips) dimensions: sets x trials x flies x conditions x taus
load([data_directory data_filename '_phistar.mat']);
phi_stars = phis;

disp('loaded')

%% Separate

phi_threes = separate_shared(phi_threes, 'phi_threes');
phi_stars = separate_shared(phi_stars, 'phi_stars');

%% Save

disp('Saving');
if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename '_phithree.mat'], 'phi_threes');
save([results_directory results_filename '_phistar.mat'], 'phi_stars');

disp('Saved');

%% Function: separate phis into sets with a sharing pair, and sets without

function [share_separated] = separate_shared(phis, phi_name)
%
%
% phis should always have the following fields:
%   nChannels, channel_sets, taus, nBins, mips, ('phi_name')
% phis should also have:
%   for phi_threes: state_counters, state_phis, tpms
%   for phi_stars: phi_stars_normalised, mis, hconds, hs, mi_stars
%
% Inputs:
%
% Outputs:
%   share_separated = phis, except with a new dimension (i.e. phis was 1x3, share_separated is 2x3);
%       the new dimension holds the phis for channel sets which have at least one pair of channels which
%       share and electrode

share_separated = cell(2, length(phis));
for nChannels_counter = 1 : length(phis)
    
    share_vector = zeros(size(phis{nChannels_counter}.channel_sets, 1), 1);
    
    % iterate through each channel set and build a logical vector (1s if set has a pair which shares and electrode)
    for set = 1 : size(phis{nChannels_counter}.channel_sets, 1)
        share_vector(set) = share_pair(phis{nChannels_counter}.channel_sets(set, :));
    end
    share_vector = logical(share_vector);
    
    % Common struct fields
    share_separated{1, nChannels_counter}.nChannels = phis{nChannels_counter}.nChannels;
    share_separated{1, nChannels_counter}.taus = phis{nChannels_counter}.taus;
    share_separated{1, nChannels_counter}.nBins = phis{nChannels_counter}.nBins;
    share_separated{2, nChannels_counter}.nChannels = phis{nChannels_counter}.nChannels;
    share_separated{2, nChannels_counter}.taus = phis{nChannels_counter}.taus;
    share_separated{2, nChannels_counter}.nBins = phis{nChannels_counter}.nBins;
    
    % Channel sets which do not have a sharing pair
    share_separated{1, nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets(not(share_vector), :);
    share_separated{1, nChannels_counter}.(phi_name) = phis{nChannels_counter}.(phi_name)(not(share_vector), :, :, :, :);
    
    % Channel sets which do have a sharing pair
    share_separated{2, nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets(share_vector, :);
    share_separated{2, nChannels_counter}.(phi_name) = phis{nChannels_counter}.(phi_name)(share_vector, :, :, :, :);
    
    % phi_threes and phi_stars have different mip structures
    if strcmp(phi_name, 'phi_threes')
        share_separated{1, nChannels_counter}.mips = phis{nChannels_counter}.mips(:, not(share_vector), :, :, :);
        share_separated{2, nChannels_counter}.mips = phis{nChannels_counter}.mips(:, share_vector, :, :, :);
    else %if strcmp(phi_name, 'phi_stars')
        share_separated{1, nChannels_counter}.mips = phis{nChannels_counter}.mips(not(share_vector), :, :, :, :);
        share_separated{2, nChannels_counter}.mips = phis{nChannels_counter}.mips(share_vector, :, :, :, :);
    end
    
    % Extra phi_three fields
    if isfield(phis{nChannels_counter}, 'state_counters')
        share_separated{1, nChannels_counter}.state_counters = phis{nChannels_counter}.state_counters(:, not(share_vector), :, :, :, :);
        share_separated{2, nChannels_counter}.state_counters = phis{nChannels_counter}.state_counters(:, share_vector, :, :, :, :);
    end
    if isfield(phis{nChannels_counter}, 'state_phis')
        share_separated{1, nChannels_counter}.state_phis = phis{nChannels_counter}.state_phis(:, not(share_vector), :, :, :);
        share_separated{2, nChannels_counter}.state_phis = phis{nChannels_counter}.state_phis(:, share_vector, :, :, :);
    end
    if isfield(phis{nChannels_counter}, 'tpms')
        share_separated{1, nChannels_counter}.tpms = phis{nChannels_counter}.tpms(:, :, not(share_vector), :, :, :);
        share_separated{2, nChannels_counter}.tpms = phis{nChannels_counter}.tpms(:, :, share_vector, :, :, :);
    end
    
    % Extra phi_star fields
    if isfield(phis{nChannels_counter}, 'phi_stars_normalised')
        share_separated{1, nChannels_counter}.phi_stars_normalised = phis{nChannels_counter}.phi_stars_normalised(not(share_vector), :, :, :, :);
        share_separated{2, nChannels_counter}.phi_stars_normalised = phis{nChannels_counter}.phi_stars_normalised(share_vector, :, :, :, :);
    end
    if isfield(phis{nChannels_counter}, 'mis')
        share_separated{1, nChannels_counter}.mis = phis{nChannels_counter}.mis(not(share_vector), :, :, :, :);
        share_separated{2, nChannels_counter}.mis = phis{nChannels_counter}.mis(share_vector, :, :, :, :);
    end
    if isfield(phis{nChannels_counter}, 'hconds')
        share_separated{1, nChannels_counter}.hconds = phis{nChannels_counter}.hconds(not(share_vector), :, :, :, :);
        share_separated{2, nChannels_counter}.hconds = phis{nChannels_counter}.hconds(share_vector, :, :, :, :);
    end
    if isfield(phis{nChannels_counter}, 'hs')
        share_separated{1, nChannels_counter}.hs = phis{nChannels_counter}.hs(not(share_vector), :, :, :, :);
        share_separated{2, nChannels_counter}.hs = phis{nChannels_counter}.hs(share_vector, :, :, :, :);
    end
    if isfield(phis{nChannels_counter}, 'mi_stars')
        share_separated{1, nChannels_counter}.mi_stars = phis{nChannels_counter}.mi_stars(not(share_vector), :, :, :, :);
        share_separated{2, nChannels_counter}.mi_stars = phis{nChannels_counter}.mi_stars(share_vector, :, :, :, :);
    end
end

end

%% Function: for a channel set, determine if there is a pair which shares an electrode

function [sharing_pair] = share_pair(channel_set)
% Given a set of channels, determines if there is at least one pair which shares an electrode
%
% Inputs:
%   channel_set = vector of channel labels (actual channel labels which correspond to electrode positions)
% Outputs:
%   sharing_pair = 1 if there is a pair of channels which shares an electrode, 0 otherwise

sharing_pair = 0;

for counter_a = 1 : length(channel_set)
    % Start 2nd iteration from counter_a to avoid researching (i.e. a,b is the same as b,a)
    for counter_b = counter_a+1 : length(channel_set)
        if abs(channel_set(counter_a) - channel_set(counter_b)) == 1
            sharing_pair = 1;
            return
        end
    end
end

end

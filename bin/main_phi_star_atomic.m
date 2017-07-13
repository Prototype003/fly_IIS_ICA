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

flies = (1:13); % This will be parallelised
nChannels = (10:15); % This determines how many channels to consider at a time
nBins = [1]; % Script currently only supports the case of 1 bin (i.e. epoch length - tau)
taus = [4, 8, 16];

nFlies = length(flies);

data_directory = 'workspace_results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

results_directory = 'results/';
results_filename = [data_file '_detrend' num2str(prep_detrend) '_zscore' num2str(prep_zscore) '_nChannels' num2str(nChannels(1)) 't' num2str(nChannels(end)) '_phistar_atomic'];

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
    
    % channels sets x trials x flies x conditions x taus x nBins
    phi_stars = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    phi_stars_normalised = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    mis = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    hconds = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    hs = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    mi_stars = zeros(size(channel_sets, 1), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5), length(taus), length(nBins));
    
    for fly_counter = 1 : nFlies
        fly = flies(fly_counter);
        disp(['nChannels' num2str(channels) ' fly' num2str(fly)]); 
        for nBins_counter = 1 : length(nBins)
            for tau_counter = 1 : length(taus)
                tau = taus(tau_counter);
                for channel_set = 1 : size(channel_sets, 1)
                    for condition = 1 : size(fly_data, 5)
                        for trial = 1 : size(fly_data, 3)
                            %disp(['fly' num2str(fly) ' tau' num2str(taus(tau_counter)), ' combo' num2str(channel_set) ' condition' num2str(condition) ' trial' num2str(trial)]); % Keep in mind this is slow
                            % Compute covariances
                            [cov_present_present, cov_present_past, cov_past_past] = Cov_comp_sample(fly_data(:, channel_sets(channel_set, :), trial, fly, condition)', tau);
                            
                            % Compute phi
                            [phi_tmp,...
                                ~, ~, ~,...
                                mis(channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                ~, ~,...
                                hconds(channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                hs(channel_set, trial, fly, condition, tau_counter, nBins_counter),...
                                mi_stars(channel_set, trial, fly, condition, tau_counter, nBins_counter)] = phi_comp(cov_present_present, cov_present_past, cov_past_past);
                            
                            phi_stars(channel_set, trial, fly, condition, tau_counter, nBins_counter) = phi_tmp(1);
                            phi_stars_normalised(channel_set, trial, fly, condition, tau_counter, nBins_counter) = phi_tmp(2);
                            
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
    
    phis{nChannels_counter}.nChannels = channels;
    phis{nChannels_counter}.channel_sets = channel_sets;
    phis{nChannels_counter}.taus = taus;
    phis{nChannels_counter}.nBins = nBins;
    
end

%% SAVE

if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename '.mat'], 'phis');
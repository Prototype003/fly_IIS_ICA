%% Description

%{
For comparing phi-three values with correlation values
%}

%% Setup

measure = 'phi_three';
global_tpm = 0;
bin_location = '';

addpath('figure_code/');

%% Load phi

% Load phi-3 values
%load('phi_3/results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM0.mat');

[phis, measure_string] = phi_load(measure, global_tpm, bin_location);

%% Load raw data

raw = load('workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat');
fly_data = raw.fly_data;

%% Calculate correlations between raw channel data

% Inefficient, because pairwise correlations are costantly recomputed
% Roughly 4 mins
corrs = cell(size(phis));
for nChannels = 1 : length(phis)
    tic;
    channel_sets = phis{nChannels}.channel_sets;
    
    phi_dims = size(phis{nChannels}.phis);
    r_means = zeros(phi_dims(1:end-1));
    r_stds = zeros(size(r_means));
    rsquared_means = zeros(size(r_means));
    rsquared_stds = zeros(size(r_means));
    
    for fly = 1 : size(fly_data, 4)
        for condition = 1 : size(fly_data, 5)
            for trial = 1 : size(fly_data, 3)
                
                for set_counter = 1 : size(channel_sets, 1)
                    channels = channel_sets(set_counter, :);
                    
                    rs = [];
                    r_counter = 1;
                    for channel_a = 1 : length(channels)
                        for channel_b = channel_a + 1 : length(channels)
                            
                            a = squeeze(fly_data(:, channels(channel_a), trial, fly, condition));
                            b = squeeze(fly_data(:, channels(channel_b), trial, fly, condition));
                            
                            r = corrcoef(a, b);
                            rs(r_counter) = r(1, 2);
                            r_counter = r_counter + 1;
                        end
                    end
                    
                    rsquareds = rs.^2;
                    
                    r_means(set_counter, trial, fly, condition) = mean(rs);
                    r_stds(set_counter, trial, fly, condition) = std(rs);
                    rsquared_means(set_counter, trial, fly, condition) = mean(rsquareds);
                    rsquared_stds(set_counter, trial, fly, condition) = std(rsquareds);
                    
                end
                
            end
        end
    end
    
    corrs{nChannels}.r_means = r_means;
    corrs{nChannels}.r_stds = r_stds;
    corrs{nChannels}.rsquared_means = rsquared_means;
    corrs{nChannels}.rsquared_stds = rsquared_stds;
    toc
end

%% Check for correlation between correlations and phi values

flies = [1 2 3 4 5 6 7 8 9 10 12 13];
figure;
phi_vals = phis{3}.phis(:, :, flies, 1, 1);
corr_vals = corrs{3}.rsquared_means(:, :, flies, 1);
scatter(phi_vals(:), corr_vals(:), '.');
xlabel('Phi'); ylabel('rsquared');
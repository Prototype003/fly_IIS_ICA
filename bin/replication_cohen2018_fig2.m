%% Description

%{

Replication of Dror's 2018 paper - Isoflurane Impairs Low-Frequency
Feedback but Leaves High-Frequency Feedforward Connectivity Intact in the
Fly Brain (eNeuro)

Replication of figure 2

%}

%% Setup

flies = (1:13);
conditions = (1:2);
channels = (1:15); % Note - Dror excludes outermost channel

% Incorrect
periph = (2:7); % "following the division of channels into peripheral channels=[1–6] and central channels=[9–14]"
center = (10:15);

% Correct
periph = (9:14);
center = (1:6);

chronux_params = struct();
chronux_params.tapers = [2 3]; % 3 tapers, 3 = 2*2 - 1
chronux_params.Fs = 1000; % 1000 Hz downsampled sampling rate of the input data
chronux_params.fpass = [0 50]; % This is the frequency range we're interested in
chronux_params.pad = 1; % I think Dror used 1, apparently higher padding gives higher frequency resolution (more frequency bins) but doesn't affect calculation
% .err - default no errorbars
% .trialave - default no averaging

nFlies = length(flies);
nConditions = length(conditions);
nChannels = length(channels);

across_flies = 1;

data_directory = 'workspace_results/';
data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat';

results_directory = 'workspace_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across' num2str(across_flies) '.mat'];
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_medianSplit_power_classification_across' num2str(across_flies) '.mat'];

%% Load

disp('loading');
load([data_directory data_filename]);
disp('loaded');

%% Preprocessing
% If required

% Detrend (no z-score)
for condition = 1 : size(fly_data, 5)
    for fly = 1 : size(fly_data, 4)
        for channel = 1 : size(fly_data, 2)
            for trial = 1 : size(fly_data, 3)
                fly_data(:, channel, trial, fly, condition) = detrend(fly_data(:, channel, trial, fly, condition));
            end
        end
    end
end


% 3 tapers, half-bandwidth 0.89 Hz



% Median split
% fly_meds = median(fly_data, 1);
% fly_meds = repmat(fly_meds, [size(fly_data, 1) 1 1 1 1]);
% fly_data = fly_data > fly_meds;

%% Calculate power
% Chronux continuous data - first dimension is time, second dimension is trials/channels

disp('Calculating power');

if across_flies >= 0
    % Preserve trials
    
    % Results matrix: frequencies x trials x channels x condition x fly
    powers = zeros(410, size(fly_data, 3), nChannels, nConditions, nFlies);
    
    for condition = 1 : nConditions
        for fly_counter = 1 : nFlies
            fly = flies(fly_counter);
            
            for channel_counter = 1 : nChannels
                channel = channels(channel_counter);
                
                % Get trials for the channel (input dimensions to chronux function should be time x trials
                trials = permute(fly_data(:, channel, :, fly, condition), [1 3 2 4 5]);
                
                % Compute power spectrums
                [spectrums, frequencies] = mtspectrumc(trials, chronux_params);
                
                % Store
                powers(:, :, channel_counter, condition, fly_counter) = spectrums;
            end
            
        end
    end
    
else % across == 1
    % Calculate across concatenated trials
    
    % Results matrix: frequencies x trials x channels x condition x fly
    powers = zeros(3277, 1, nChannels, nConditions, nFlies);
    
    for condition = 1 : nConditions
        for fly_counter = 1 : nFlies
            fly = flies(fly_counter);
            
            for channel_counter = 1 : nChannels
                channel = channels(channel_counter);
                
                % Get trials for the channel (input dimensions to chronux function should be time x trials
                trials = permute(fly_data(:, channel, :, fly, condition), [1 3 2 4 5]);
                
                % Concatenate trials into one big trial
                trial = zeros(size(trials, 1) * size(trials, 2), size(trials, 3));
                row_counter = (1:size(trials, 1));
                for trial_counter = 1 : size(trials, 2)
                    trial(row_counter, :) = trials(:, trial_counter, :);
                    row_counter = row_counter + size(trials, 1);
                end
                
                % Compute power spectrums
                [spectrums, frequencies] = mtspectrumc(trial, chronux_params);
                
                % Store
                powers(:, :, channel_counter, condition, fly_counter) = spectrums;
            end
            
        end
    end
    
end

powers = 10*log10(powers);

disp('log powers calculated');

%% Make figure

condition = 2;

figure;

% Mean across trials, channels, flies
plot(frequencies, squeeze(mean(mean(mean(powers(:, :, center, condition, :), 2), 3), 5)));
hold on;
plot(frequencies, squeeze(mean(mean(mean(powers(:, :, periph, condition, :), 2), 3), 5)));
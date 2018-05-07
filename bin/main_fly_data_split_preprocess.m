%% DESCRIPTION

%{

Carries out pre-processing of fly data in preparation for highly parallel phi-3 computation

Pre-processing: each trial is binarised based on median split

Data is split into - trial, fly, condition

%}

%% Setup

channels = (1:15);
nChannels = 4;

channel_sets = nchoosek(channels, nChannels) - 1;
    
data_directory = 'workspace_results/';
data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat';

results_directory = 'fly_data_split/';
results_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

if ~isdir(results_directory)
    mkdir(results_directory)
end

%% Load

load([data_directory data_filename]);

%% Pre-process

% Dimensions of fly data
[nSamples, ~, nTrials, nFlies, nConditions] = size(fly_data);

for fly = 1 : nFlies
    for condition = 1 : nConditions
        for trial = 1 : nTrials
            
            % Get data
            data = fly_data(:, :, trial, fly, condition);
            
            % Binarise each channel
            medians = median(data, 1);
            medians = repmat(medians, [size(data, 1) 1]);
            data = data > medians;
            
            % Save data
            results_suffix = ['_f' num2str(fly) 'c' num2str(condition) 't' num2str(trial)];
            save([results_directory results_prefix results_suffix '.mat'], 'data', 'nChannels', 'channel_sets');
            
            disp(results_suffix);
        end
    end
end
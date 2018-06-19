%% DESCRIPTION

%{

Takes split fly data and gives a new channel_sets matrix

%}

%% Setup

channels = (1:15);
nChannels_new = 3;

channel_sets_new = nchoosek(channels, nChannels_new) - 1; % Python indexing

source_directory = 'fly_data_split/';
source_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

dest_directory = 'fly_data_split_3ch/';
dest_prefix = source_prefix;

if ~isdir(dest_directory)
    mkdir(dest_directory)
end

%% Load and update each data file

% Data dimensions
nFlies = 13;
nConditions = 2;
nTrials = 8;

for fly = 1 : nFlies
    for condition = 1 : nConditions
        for trial = 1 : nTrials
            
            file_id = ['_f' num2str(fly) 'c' num2str(condition) 't' num2str(trial)];
            
            % Load data file
            load([source_directory source_prefix file_id '.mat']);
            
            % Update
            channel_sets = channel_sets_new;
            nChannels = nChannels_new;
            
            % Save into new location/file
            save([dest_directory source_prefix file_id '.mat'], 'data', 'nChannels', 'channel_sets');
            
        end
    end
end

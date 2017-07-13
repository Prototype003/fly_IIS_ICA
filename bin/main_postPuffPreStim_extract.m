%% DESCRIPTION

%{
This takes takes the Dror's data at the step where line noise has been
removed, and extracts the data between the puff and the stimulus block.

Source data: /export/kani/shared/drorcNew/Current/Data/Flies/Processed/trials_anesthDescAdded23092014_endogTrialsIncRecv_pffDcmp2_endTrn_Amnd_split2250_bPlrRerefTyp1_lnNsRmvd

The data at this stage has been epoched and bipolar rereferenced

Reformats data into the following form:
    Matrix - sample(2250) x channel(15) x trial(8) x fly(13) x condition(2)
%}

%% SETUP

data_directory = '../../flies/fly_data_lineNoiseRemoved/';
data_file = 'trials.mat';

results_directory = 'workspace_results/';
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

flies = (1:13);
nFlies = length(flies);

%% Get postpuff-prestim blocks

fly_ids = cell(nFlies, 1);
relevant_periods = cell(nFlies, 1);
for fly_counter = 1 : nFlies
    fly = flies(fly_counter);
    fly_id = fly_names(fly);
    load([data_directory fly_id '/' data_file]);
    
    fly_ids{fly_counter} = fly_id;
    
    % 'Fix' for one fly which had a longer postpuff-prestim period,
    % resulting in more than 8 2.25s epochs - take only the last 8
    if strcmp(fly_id, 'Analyzed_WalkingDrosDror145707072014')
        trials(50:57) = [];
    end
    
    % Extract trials/periods with an air puff
    relevant_periods{fly_counter} = extract_spontaneous_data(trials);
    
end

trials = struct('fly_id', fly_ids, 'relevant_data', relevant_periods);

% Reformat to matrix
fly_data = matrix_format(trials);

%% Save

if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename '.mat'], 'fly_data');


%% FUNCTIONS

function [relevant_trials] = extract_spontaneous_data(trials)
% Extracts entries from the non-scalar struct array trials where
% puffStartStop specifies a postpuff-prestim block
%
% Inputs:
%   trials = non-scalar struct array, must have field puffStartStop
%
% Outputs:
%   puff_trials = non-scalar struct array

relevant_trial_counter = 1;
for trial = 1 : length(trials)
    if strcmp(trials(trial).anesthNm, '0_postPuffPreStim') || strcmp(trials(trial).anesthNm, '2_postPuffPreStim')
        relevant_trials(relevant_trial_counter) = trials(trial);
        relevant_trial_counter = relevant_trial_counter + 1;
    end
end

end

function [reformatted] = matrix_format(trials)
% Reformats the nonscalar struct array trials into a matrix with dimensions:
% sample(2250) x channel(15) x trial(8) x fly(13) x condition(2)
% Condition 1 = air; condition 2 = anaesthesia
%
% Assumes 2 conditions, and equal number of trials for each condition (same
% number of trials per fly, per condition) - due to matrix restrictions,
% dimensions have to be the same for all flies
%
% Inputs:
%   trials = nonscalar struct array with each entry corresponding to a fly
%   and holding its data in the field relevant_data
%
% Outputs:
%   reformatted = matrix with dimensions as described above

[nChannels, nSamples] = size(trials(1).relevant_data(1).LFP);
nTrials = length(trials(1).relevant_data) / 2;
nFlies = length(trials);

reformatted = zeros(nSamples, nChannels, nTrials, nFlies, 2);

for fly = 1 : nFlies
    condition_counters = [0 0];
    for trial = 1 : length(trials(1).relevant_data)
        
        if strcmp(trials(fly).relevant_data(trial).anesthNm, '0_postPuffPreStim')
            condition = 1;
        elseif strcmp(trials(fly).relevant_data(trial).anesthNm, '2_postPuffPreStim')
            condition = 2;
        else
            disp('ERROR: More than 2 conditions');
        end
        
        condition_counters(condition) = condition_counters(condition) + 1;
        reformatted(:, :, condition_counters(condition), fly, condition) = trials(fly).relevant_data(trial).LFP';
        
    end
end

end

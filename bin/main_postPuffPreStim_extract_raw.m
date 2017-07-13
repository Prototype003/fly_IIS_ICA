%% DESCRIPTION

%{
This script takes the relevant extracts data corresponding to rest periods
(i.e. periods which have an air-puff) - essentially takes the whole dataset
and saves only the data corresponding to the relevant timeblocks

It saves two outputs:
    Periods between an air puff and a stimulus block (in that order)
    Periods between a stimulus block and an air puff (in that order)

Outputs are not processed (the data is still raw)
%}

%% SETUP

data_directory = '../../flies/fly_data/half_brain_flies/';
data_file = 'trials.mat';

results_directory = 'workspace_results/';
results_filename = 'postPuffpreStim';

flies = (1:13);
nFlies = length(flies);

%% Get air puff periods

fly_ids = cell(nFlies, 1);
puff_periods = cell(nFlies, 1);
for fly_counter = 1 : nFlies
    fly = flies(fly_counter);
    fly_id = fly_names(fly);
    load([data_directory fly_id '/' data_file]);
    
    % We need to apply manual fixes (as identified in grabEndogTrials.m) to
    % two flies
    if strcmp(fly_id, 'Analyzed_WalkingDrosDror145707072014')
        concatLFP=[trials(321).LFP trials(322).LFP trials(323).LFP];
        trials(321).LFP=concatLFP;
        trials(321).endTime=trials(323).endTime;
        trials(321).duration=size(concatLFP,2);
        %trials(321).endTime-trials(321).startTime % This seems to be a check printed out for the user (total duration)
    elseif strcmp(fly_id, 'Analyzed_WalkingDrosDror160014072014')
        concatLFP=[trials(161).LFP trials(162).LFP trials(163).LFP];
        trials(161).LFP=concatLFP;
        trials(161).endTime=trials(163).endTime;
        trials(161).duration=size(concatLFP,2);
        %trials(161).endTime-trials(161).startTime; % This seems to be a check printed out for the user (total duration)
        trials(161).puffStartStop(2,:)= trials(163).puffStartStop(1,:);
    end
    
    fly_ids{fly_counter} = fly_id;
    
    % Extract trials/periods with an air puff
    puff_periods{fly_counter} = extract_spontaneous_data(trials);
    
    % Extract only pre-stimulus puff trials (for trials with multiple puffs)
    puff_periods{fly_counter} = extract_prestim_puffs(puff_periods{fly_counter});
    
    % Extract samples corresponding to the period between the puff and the
    % stimulus block
    puff_periods{fly_counter} = extract_prestim_samples(puff_periods{fly_counter});
    
end

trials = struct('fly_id', fly_ids, 'puff_data', puff_periods);

%% Save

if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename '.mat'], trials);


%% FUNCTIONS

function [puff_trials] = extract_spontaneous_data(trials)
% Extracts entries from the non-scalar struct array trials where
% puffStartStop is not NaN
%
% Inputs:
%   trials = non-scalar struct array, must have field puffStartStop
%
% Outputs:
%   puff_trials = non-scalar struct array

puff_trial_counter = 1;
for trial = 1 : length(trials)
    if (length(trials(trial).puffStartStop) > 1) &&... % If length is one, then there is a NaN (puff periods have a puff start and stop, so length is at least 2)
            (trials(trial).duration > 2500) &&... % To exclude air puffs during stimulus blocks (in which case trials are ~2400ms)
            (...
            strcmp(trials(trial).anesthNm, '0') ||... % air
            strcmp(trials(trial).anesthNm, '0trn2') ||... % air to full anaesthesia
            strcmp(trials(trial).anesthNm, '1.5trn2') ||... % graded anaesthesia to full anaesthesia (should be equivalent to 'air to full')
            strcmp(trials(trial).anesthNm, '2trn0')... % full anaesthesia to air (rest)
            )
        puff_trials(puff_trial_counter) = trials(trial);
        puff_trial_counter = puff_trial_counter + 1;
    end
end

end

function [puff_data] = extract_prestim_puffs(puff_data)
% Removes all puffs except for the last one (the one corresponding to
% postpuff-prestim
% Assumes that the chronologically last puff is in the last row of
% puffStartStop
%
% Inputs:
%   puff_data = non-scalar struct array with field puffStartStop
%   (puffStartStop is a 2D matrix, with each row corresponding to a puff)
%
% Outputs:
%   puff_data = inputted puff_data, but puffStartStop should all be 1D,
%   corresponding to a single puff

for trial = 1 : length(puff_data)
    puff_data(trial).puffStartStop = puff_data(trial).puffStartStop(end, :);
end

end

function [puff_data] = extract_prestim_samples(puff_data)
% Keeps samples which are are between the puff stop and the end of the
% sampling block (i.e. between the last of puffStartStop and endTime)
%
% Inputs:
%   puff_data = non-scalar struct arrawy with field puffStartStop, endTime,
%   and LFP. The last element of puffStartStop should be the end of the
%   puff
%
% Outputs:
%   puff_data = inputted puff_data, but LFP samples outside of the
%   postpuff-prestim period are removed

for trial = 1 : length(puff_data)
    block_duration = puff_data(trial).endTime - (puff_data(trial).puffStartStop(end)); % Duration between puff and stimulus (-1)
    
    % Restrict LFP time range
    puff_data(trial).LFP = puff_data(trial).LFP(:, size(puff_data(trial).LFP, 2)-block_duration : end);
    % Update startTime to match new LFP range (endTime should be the same)
    puff_data(trial).startTime = puff_data(trial).endTime - block_duration;
end

% Duration is not needed (and a little misleading, as it is the number
% of samples - 1 (e.g. startTime - endTime : 100 - 1 = 99, but duration
% is 100, not 99)
puff_data = rmfield(puff_data, 'duration');

end

%%
% Splitting into 2.25s epochs - start from the end of the 18s period (not
% from the air puff itself - justification is that this will prioritise
% discarding data with residual air flow as a potential stimulant)

% Furthermore we will need to try and do this and that as well. The reasons
% for this are not so interesting, so we won't go through them. Further
% motivation for this as inspired by Dror is to avoid forcing certain
% predetermined outcomes



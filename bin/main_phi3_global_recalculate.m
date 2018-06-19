%% Description

%{

This is for recalculating phi-3 using the full 18s TPM

Output - the trial dimension is collapsed to 1 trial, with its phi value
being the average phi across states, after weighting for occurrence of each
state across all trials

Output dimensions will be the same - only the phi-value matrix is updated to
only have one trial

%}

%% Load

source_directory = 'results/';
source_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree.mat';

dest_directory = 'results/';
dest_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree_global.mat';

disp('loading');
load([source_directory source_file]);
disp('loaded');

%% Recalculate based on state counters

for nChannels = 1 : length(phis)
    
    % Sum state_counters across trials
    state_counters_global = permute(sum(phis{nChannels}.state_counters, 3), [1 2 4 5 6 3]);
    
    % Multiply phis for each state by the number of times they occurred
    summed_phis = phis{nChannels}.state_phis .* state_counters_global;
    
    % Sum multiplied phis across states
    summed_phis = permute(sum(summed_phis, 1), [2 3 4 5 1]);
    
    % Divide summed phis by total number of states
    phis{nChannels}.phi_threes = permute(summed_phis ./ permute(sum(state_counters_global, 1), [2 3 4 5 1]), [1 5 2 3 4]); % And add in the singleton trial dimension
end

%% Save

save([dest_directory dest_file], 'phis');

disp('saved');
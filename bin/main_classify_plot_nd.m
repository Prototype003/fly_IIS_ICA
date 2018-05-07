%% Description

%{
Plot values as a function of channel location.

Plots are 2D for 2-channels, 3D for 3-channels, and 4D for 4-channels
%}

%% Settings

nChannels = 2;

results_directory = 'workspace_results/';
nConditions = 2;

accuracies
across = 1;

if across == 1
    results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
else
    results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
end

%% Load

load([results_directory results_filename]);

%% 2 channel plot

figure;

values = accuracies{
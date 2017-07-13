%% DESCRIPTION

%{
This script tests the first hypothesis: phi in the air condition is greater than phi in the isoflurane condition

Tests for all nChannels in loaded datafile

%}

%% SETUP

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_phithree'...
    ];

%% LOAD

load([data_directory data_filename '.mat']);

%% t-test (after averaging across tau, channel sets, and trials)

nChannels = length(phis);

averaged_results = cell(nChannels, 1);

for nChannels_counter = 1 : nChannels
    phi_struct = phis{nChannels_counter};
    channels_used = phi_struct.nChannels;
    averaged_results{nChannels_counter}.nChannels = channels_used;
    
    % Does it matter what order we average in?
    % i.e. across tau, then channel sets, then trials vs across channel sets, tau, then trials?
    phi_values = phi_struct.phi_threes;
    phi_values = squeeze(mean(mean(mean(phi_values, 5), 1), 2));
    
    [averaged_results{nChannels_counter}.H,...
        averaged_results{nChannels_counter}.P,...
        averaged_results{nChannels_counter}.CI,...
        averaged_results{nChannels_counter}.stats] = ttest(phi_values(:, 1), phi_values(:, 2));
end

%% t-test without averaging (at each tau)

%% Taking into account nesting

nChannels = length(phis);

nested_results = cell(nChannels, 1);
for nChannels_counter = 1 : nChannels
    phi_struct = phis{nChannels_counter};
    channels_used = phi_struct.nChannels;
    nested_results{nChannels_counter}.nChannels = channels_used;
    
    % We are using a one-way ANOVA : main effect of condition
    % Or are we using a two-way ANOVA : main effects of condition and tau? (or can we treat tau as nested?)
    condition_grouping = [];
    fly_grouping = [];
    trial_grouping = [];
    tau_grouping = [];
    % Nesting order:
    % fly > trial > channel combination > tau
    
    % Channel combinations are nested within trials, which are nested within flies
    
    
    % Reformat phi matrix into a vector while allowing for grouping/nesting
    phi_vector = zeros(numel(phi_struct.phi_threes), 1);
    
    
end

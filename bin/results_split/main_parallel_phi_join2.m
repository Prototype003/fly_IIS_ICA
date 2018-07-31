%% DESCRIPTION

%{

For joining highly parallelised phi3 results (4 channels) into one file per fly

%}

%% Setup

nFlies = 13;
nConditions = 2;
taus = [4];
nSets = 1365;
nTrials = 8;
channels = (1:15);
nChannels = 4;
nStates = nChannels^2;
nPartitions = 14; % Number of bipartitions for a set of 4

data_detrended = 0;
data_zscored = 0;

data_directory_prefix = ''; % Script goes directly in the data directory
data_file_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_phithree_allPartitions_';
data_file_suffix = '.mat';

results_directory = ''; % Results are saved directly in the data directory

%% Concatenated results structure

phis = cell(1);

results = struct();
results.nChannels = nChannels;
results.channel_sets = nchoosek(channels, nChannels);
results.taus = taus;

results.phi_threes = zeros(nSets, nTrials, nConditions, length(taus));
results.state_counters = zeros(nStates, nSets, nTrials, nConditions, length(taus));
results.state_phis = zeros(nStates, nSets, nTrials, nConditions, length(taus));
results.tpms = zeros(nStates, nStates, nSets, nTrials, nConditions, length(taus));
results.mips = cell(nStates, nSets, nTrials, nConditions, length(taus));
results.state_partitions = cell(nStates, nPartitions, nSets, nTrials, nConditions, length(taus));
results.state_partitions_phis = zeros(nStates, nPartitions, nSets, nTrials, nConditions, length(taus));

%% Concatenate results
for fly = 1 : nFlies
    disp(['fly' num2str(fly)]);
    data_directory = [data_directory_prefix 'fly' num2str(fly, '%02.f') '/'];
    for condition = 1 : nConditions
        for tau_counter = 1 : length(taus)
            disp(['fly' num2str(fly) '_cond' num2str(condition) '_tau' num2str(taus(tau_counter))]);
            for set = 1 : nSets
                for trial = 1 : nTrials
                    data_file_infix = [...
                        'f' num2str(fly, '%02.f')...
                        'c' num2str(condition)...
                        'tau' num2str(taus(tau_counter))...
                        's' num2str(set-1, '%04.f')... % set-1 because of python indexing/naming
                        't' num2str(trial)...
                        ];
                    data_file = [data_file_prefix data_file_infix data_file_suffix];
                    
                    load([data_directory data_file]);
                    
                    results.phi_threes(set, trial, condition, tau_counter) = phi.phi;
                    results.state_counters(:, set, trial, condition, tau_counter) = phi.state_counters;
                    results.state_phis(:, set, trial, condition, tau_counter) = phi.state_phis;
                    results.tpms(:, :, set, trial, condition, tau_counter) = phi.tpm;
                    results.mips(:, set, trial, condition, tau_counter) = phi.mips;
                    results.state_partitions(:, :, set, trial, condition, tau_counter) = phi.state_partitions;
                    results.state_partitions_phis(:, :, set, trial, condition, tau_counter) = phi.state_partitions_phis;
                end
            end
        end
    end
    
    phis{1} = results;
    
    results_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels4t' num2str(nChannels)...
    '_phithree'...
    '_nonGlobal'...
    '_allPartitions'...
    'fly' num2str(fly)...
    '.mat'...
    ];
    
    save([results_directory results_file], 'phis', '-v7.3');
end

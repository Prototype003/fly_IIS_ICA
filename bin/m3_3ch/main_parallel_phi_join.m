%% DESCRIPTION

%{

For joining highly parallelised phi3 results (4 channels) into one file (all flies)

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
results_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels4t' num2str(nChannels)...
    '_phithree'...
    '_nonGlobal'...
    '.mat'
    ];

%% Concatenated results structure

phis = cell(1);

results = struct();
results.nChannels = nChannels;
results.channel_sets = int8(nchoosek(channels, nChannels));
results.taus = taus;

phi_threes = zeros(nSets, nTrials, nFlies, nConditions, length(taus));
state_counters = zeros(nStates, nSets, nTrials, nFlies, nConditions, length(taus));
state_phis = zeros(nStates, nSets, nTrials, nFlies, nConditions, length(taus));
tpms = zeros(nStates, nStates, nSets, nTrials, nFlies, nConditions, length(taus));
mips = cell(nStates, nSets, nTrials, nFlies, nConditions, length(taus));
%state_partitions = cell(nStates, nPartitions, nSets, nTrials, nFlies, nConditions, length(taus));
%state_partitions_phis = zeros(nStates, nPartitions, nSets, nTrials, nFlies, nConditions, length(taus));

%% Concatenate results
for fly = 1 : nFlies
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
                    
                    phi_threes(set, trial, fly, condition, tau_counter) = single(phi.phi);
                    state_counters(:, set, trial, fly, condition, tau_counter) = int16(phi.state_counters);
                    state_phis(:, set, trial, fly, condition, tau_counter) = single(phi.state_phis);
                    tpms(:, :, set, trial, fly, condition, tau_counter) = single(phi.tpm);
                    mips(:, set, trial, fly, condition, tau_counter) = mips_to_int(phi.mips);
                    %state_partitions(:, :, set, trial, fly, condition, tau_counter) = partitions_to_int(phi.state_partitions);
                    %state_partitions_phis(:, :, set, trial, fly, condition, tau_counter) = single(phi.state_partitions_phis);
                end
            end
        end
    end
end

phis{1} = results;

%% Save

whos

disp('Saving');
save([results_directory results_file])%, '-v7.3', '-nocompression');
disp('Saved');

%% Conversions for MIPS and partitions

function [mips] = mips_to_int(mips)

for state = 1 : length(mips)
    state_mip = mips{state};
    if iscell(state_mip)
        for part = 1 : length(state_mip)
            state_mip{part} = int8(state_mip{part});
        end
    else
        state_mip = int8(state_mip);
    end
    mips{state} = state_mip;
end

end

function [state_partitions] = partitions_to_int(state_partitions)

nStates = size(state_partitions, 1);
nPartitions = size(state_partitions, 2);

for state = 1 : nStates
    for partition = 1 : nPartitions
        state_partition = state_partitions{state, partition};
        if iscell(state_partition)
            for part = 1 : length(state_partition)
                state_partition{part} = int8(state_partition{part});
            end
        else
            state_partition = int8(state_partition);
        end
        state_partitions{state, partition} = state_partition;
    end
end
end
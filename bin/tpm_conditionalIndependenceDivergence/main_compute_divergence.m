%% Description

%{

Computes state-by-state TPMs, converts to state-by-node, and converts back

Finds divergence between original and back-converted SBS TPMs

%}

%% Setup

nChannels = 4;
nValues = 2;
tau = 4;

addpath('../');

% Load fly data
tic;
load('../workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat');
toc

%% Compute TPMs
% ~10 minutes on 8 cores

networks = nchoosek((1:size(fly_data, 2)), nChannels);

data = binarise_median(fly_data);

dims = size(fly_data);
tpms_orig = zeros([nValues^nChannels nValues^nChannels size(networks, 1) dims(3:end)]);
tpms_sbn = zeros([nValues^nChannels nChannels size(networks, 1) dims(3:end)]);
tpms_back = zeros(size(tpms_orig));
tpm_divergence = zeros([nValues^nChannels dims(3:end)]);
for condition = 1 : size(fly_data, 5)
    for fly = 1 : size(fly_data, 4)
        for trial = 1 : size(fly_data, 3)
            tic;
            parfor network = 1 : size(networks, 1)
                
                tpm_sbs = build_tpm(data(:, networks(network, :), trial, fly, condition), tau, nValues);
                tpms_orig(:, :, network, trial, fly, condition) = tpm_sbs;
                
                tpm_sbn = tpm_sbs2sbn(tpm_sbs);
                tpms_sbn(:, :, network, trial, fly, condition) = tpm_sbn;
                
                tpm_back = tpm_sbn2sbs(tpm_sbn);
                tpms_back = tpm_back;
                
                divergence = corr(tpm_sbs(:), tpm_back(:));
                tpm_divergence(network, trial, fly, condition) = divergence;
                
            end
            toc
        end
    end
end

%% Save

tic;
save('tpm_divergence.mat', 'tpms_orig', 'tpms_sbn', 'tpms_back', 'tpm_divergence', 'nChannels', 'networks', 'nValues', 'tau');
toc
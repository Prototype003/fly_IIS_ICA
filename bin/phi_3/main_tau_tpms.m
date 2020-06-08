%% Description

%{

Builds TPM at each tau

Then, correlate TPMs for each pair of taus

%}

%% Setup

nComponents = 4;
taus = (1:2000);

addpath('../');

filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials'...
    '_nComponents' num2str(nComponents)];

% Load ICs
load(['../ICA/results/' filename '_all.mat']);
load(['../ICA/results/' filename  '.mat']);

fly_data_original = fly_data;

%% Reformat fly data

fly_data = fly_data_original;

fly_data = permute(fly_data, [1 3 2 4 5]); % samples x trials x conditions x channels x flies

dims = size(fly_data);

% Collapse trials
fly_data = reshape(fly_data, [dims(1)*dims(2) dims(3) dims(4) dims(5)]);

%% Build TPMs from ICs
% Downsamples with offsets

binarise_method = 'median';

% Assumes (samples x channels x flies x conditions)
networks = nchoosek((1:size(fly_data, 2)), nComponents);
tpms = zeros(nComponents^2, nComponents^2, size(networks, 1), size(fly_data, 3), size(fly_data, 4), length(taus));

for tau_counter = 1 : length(taus)
    tau = taus(tau_counter);
    disp(tau);
    tic;
    for fly = 1 : size(fly_data, 3)
        for condition = 1 : size(fly_data, 4)
            for network = 1 : size(networks, 1)
                data = fly_data(:, networks(network, :), fly, condition);
                
                % Build TPM
                tmp = build_tpm_binOffsets(data, tau, 2, binarise_method);
                tpms(:, :, network, fly, condition, tau_counter) = tmp;
            end
        end
    end
    toc
end

%% Correlation TPMs between taus

% taus x taus x networks x flies x conditions
tpm_tau_corrs = zeros(size(tpms, 6), size(tpms, 6), size(tpms, 3), size(tpms, 4));

for fly = 1 : size(tpms, 4)
    tic;
    for condition = 1 : size(tpms, 5)
        for network = 1 : size(tpms, 3)
            data = reshape(tpms(:, :, network, fly, condition, :),...
                [size(tpms, 1)*size(tpms, 1) size(tpms, 6)]);
            tpm_tau_corrs(:, :, network, fly, condition) = corr(data);
        end
    end
    toc
end

%% Plot

condition = 2;

clim = [min(tpm_tau_corrs(:)) max(tpm_tau_corrs(:))];

% Plot average across flies
figure; imagesc(mean(tpm_tau_corrs(1:500, 1:500, 1, :, condition), 4), clim); colorbar;

% Plot for each fly
figure;
for fly = 1 : size(tpm_tau_corrs, 4)
    subplot(4, 4, fly);
    imagesc(tpm_tau_corrs(1:500, 1:500, 1, fly, condition), clim); colorbar;
end

%% Save

save(['tpms/' filename '_tpmTauCorr'], 'tpms', 'taus', 'tpm_tau_corrs');

%% Stats

network = 1;

sig_mat = zeros(length(taus), length(taus), size(tpm_tau_corrs, 3));
p_mat = zeros(size(sig_mat));

for t1 = 1 : size(tpm_tau_corrs, 1)
    tic;
    for t2 = t1 : size(tpm_tau_corrs, 1)
        a = squeeze(fisher_rz(tpm_tau_corrs(t1, t2, network, :, 1)));
        b = squeeze(fisher_rz(tpm_tau_corrs(t1, t2, network, :, 2)));
        
        [sig_mat(t1, t2), p_mat(t1, t2)] = ttest(a, b);
    end
    toc
end

%% Save

save(['tpms/' filename '_tpmTauCorr'], 'tpms', 'taus', 'tpm_tau_corrs', 'sig_mat', 'p_mat');

%% Plot stats

% FDR correction
q = 0.1;
p_reshape = reshape(p_mat, [numel(p_mat) 1]);
[pID pN] = FDR(p_reshape, q);
sig_fdr = reshape(p_mat < pID, size(p_mat));




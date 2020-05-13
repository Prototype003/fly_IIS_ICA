%% Description

%{
Finds correlation between TPM divergence and difference in phi
%}

%% Load divergences

load('results/tpm_divergence.mat');

%% Load phis

load('../phi_3/results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM0.mat');

%% All divergences, vs all phis

figure;
scatter(tpm_divergence(:), phis{1}.phis(:), '.');

%% Mean divergence across conditions, vs phi difference
figure;
a = mean(tpm_divergence, 4);
b = phis{1}.phis(:, :, :, 1) - phis{1}.phis(:, :, :, 2);
scatter(a(:), b(:), '.');
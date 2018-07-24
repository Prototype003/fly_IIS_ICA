%% Description

%{

Script for comparing chained phi results (chained from mean 2-ch phi)

%}

%% Common setup

results_dir = 'results_phiChain/';

%% 1ms tau vs 4ms tau

% Load
tmp1 = load([results_dir 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels5_globalTPM0_chainedPhi_tau1.mat']);
tmp2 = load([results_dir 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels5_globalTPM0_chainedPhi_tau4.mat']);
tau_1 = tmp1.phi_joined.phis;
tau_2 = tmp2.phi_joined.phis;

% Dimensions: trials x flies x conditions

% Compare all values to each other (values are matched)
[h, p] = ttest(tau_1(:), tau_2(:))

% Compare awake values to each other
a = tau_1(:, :, 1);
b = tau_2(:, :, 1);
[h, p] = ttest(a(:), b(:))

% Compare anest values to each other
a = tau_1(:, :, 2);
b = tau_2(:, :, 2);
[h, p] = ttest(a(:), b(:))

% Compare classification accuracies
tau_1 = tmp1.phi_joined.accuracies;
tau_2 = tmp2.phi_joined.accuracies;

[h, p] = ttest(tau_1, tau_2)

%% Unaveraged 4ms tau vs averaged 4ms tau

% Load
unaveraged = load([results_dir 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels5_globalTPM0_chainedPhi_tau4.mat']);
averaged = load([results_dir 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels5_globalTPM0_chainedPhi_tauBin4.mat']);
unaveraged = unaveraged.phi_joined.phis;
averaged = averaged.phi_joined.phis;

% Dimensions: trials x flies x conditions

% Compare all values to each other (values are matched)
[h, p] = ttest(unaveraged(:), averaged(:))

% Compare awake values to each other
a = unaveraged(:, :, 1);
b = averaged(:, :, 1);
[h, p] = ttest(a(:), b(:))

% Compare anest values to each other
a = unaveraged(:, :, 2);
b = averaged(:, :, 2);
[h, p] = ttest(a(:), b(:))

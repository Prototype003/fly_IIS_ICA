%% Description

%{

Compares classification results between phi and power

%}

%% Setup

results_dir = 'results/';

%% Within fly classification comparison

phi_class = load([results_dir 'phi3_svm_within.mat']);
power_class = load([results_dir 'power_svm_within.mat']);

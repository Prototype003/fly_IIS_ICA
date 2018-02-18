avi%% DESCRIPTION

%{

Paper figure 1

Example MIP search - show values for each partition

%}

%% SETUP

fly = 1;
channel_set = 1365;
condition = 1;
tau = 1;

data_directory = 'results/preformatted_results/';

% This has all set sizes (2t4), but only for one fly
data_file_star = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_medianSplit0_phistar_allPartitionsfly' num2str(fly) '.mat'];

% This has only one set_size (3), for all flies
data_file_three = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_allPartitions.mat';

%% LOAD

disp('loading');
load([data_directory data_file_star]);
phi_stars = phis;
load([data_directory data_file_three]);
phi_threes = phis;
disp('loaded');

%% Find channel set which behaves as predicted (decrease in phi during iso)

phi_diff = phi_threes{1}.phi_threes(:, :, fly, 1, 1) - phi_threes{1}.phi_threes(:, :, fly, 2, 1);
%phi_diff = phi_stars{2}.phi_stars(:, :, 1, 1, 1) - phi_stars{2}.phi_stars(:, :, 1, 2, 1);

[value, position] = max(phi_diff(:));

[max_set, trial] = ind2sub(size(phi_diff), position);

%% Plot for phi-3

% Pick a state
state = 1;

% Get partitions
partitions = phi_threes{1}.state_partitions(state, :, max_set, fly, condition, tau);

% Get phi values
partitions_phis = phi_threes{1}.state_partitions_phis(state, :, max_set, fly, condition, tau);

% Sort partitions and phi values
partition_order = partition_sort(partitions);
partitions = partitions(partition_order);
partitions_phis = partitions_phis(partition_order);

partition_labels = {...
    'A--|-->B,C',...
    'B,C--|-->A',...
    'A,B--|-->C',...
    'C--|-->A,B',...
    'A,C--|-->B',...
    'B--|-->A,C'};

figure;
subplot(1, 2, 1)
bar(partitions_phis);
set(gca, 'XTick', (1:6), 'XTickLabel', partition_labels, 'XTickLabelRotation', 45, 'YTick', [0 0.05 0.1], 'FontSize', 12);
axis([0 7 0 0.12]);
xlim([0 7]);
ylabel('\Phi', 'FontSize', 15, 'rotation', 90);

%% Plot for phi-star
figure;
% Get partitions
partitions = phi_stars{2}.partitions(:, max_set, trial, 1, condition, tau);

% Get values
partitions_phis = phi_stars{2}.partitions_mi_stars(:, max_set, trial, 1, condition, tau);

partition_labels = {...
    'AB-|-C',...
    'AC-|-B',...
    'A-|-BC',...
    'A-|-B-|-C'};

subplot(1, 2, 2);
bar(partitions_phis);
set(gca, 'XTick', (1:4), 'XTickLabel', partition_labels, 'XTickLabelRotation', 45, 'YTick', [0 0.05 0.1], 'FontSize', 12);

%axis([0 5 0 0.12]);
xlim([0 5]);
ylabel('\Phi*', 'FontSize', 15, 'rotation', 90);
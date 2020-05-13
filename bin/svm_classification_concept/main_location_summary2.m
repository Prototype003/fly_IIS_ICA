%% Description

%{

Plots and stats for seeing if there is a correlation between some set
attribute (e.g. mean channel location) and classification accuracy

%}

%% Setup

const_types = {'unpart', 'part', 'both'};
const_types = {'unpart'};
class_type = 'across';

%% Load accuracies

c_accuracies = struct(); % accuracies for each individual concept

results_dir = 'results/';
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    results_file = ['4ch_phi3Concept_' const_type '_svm_' class_type '.mat'];
    acc_tmp = load([results_dir results_file]);
    
    % Take the values at a specific cost
    cost_param = 1;
    c_accuracies.(const_type) = permute(mean(acc_tmp.cost_accuracies(:, :, :, cost_param), 2), [1 3 2]);
    
    % Take the maximum across costs
    c_accuracies.(const_type) = permute(mean(max(acc_tmp.cost_accuracies, [], 4), 2), [1 3 2]);
    
    cost_acc_all = acc_tmp.cost_accuracies;
    
end

%% Processing

% Average across all concepts within the same order

o_accuracies = c_accuracies; % accuracies averaged across concepts with the same order

max_order = acc_tmp.nChannels; % from the last loaded file, but all should have the same value

for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    
    % (networks x orders)
    o_accuracies.(const_type) = zeros(size(c_accuracies.(const_type), 1), max_order);
    
    concept_ind = 1; % 0 because we will add counters to it
    for order_c = 1 : max_order % for each order
        
        nConcepts = nchoosek(max_order, order_c);
        concept_indices = (concept_ind : concept_ind+nConcepts-1);
        
        o_accuracies.(const_type)(:, order_c) = mean(c_accuracies.(const_type)(:, concept_indices), 2);
        
        o_accuracies.(const_type)(:, order_c) = max(c_accuracies.(const_type)(:, concept_indices), [], 2);
        
        concept_ind = concept_ind + nConcepts;
        
    end
    
end

% Big phi
p_accuracies = c_accuracies.unpart(:, end); % big-phi is same for all constellation types

%% Add full constellation results

% Load classification results from using full small phi constellation

results_dir = '../svm_classification_composition/results/';

comp_accuracies = struct(); % accuracies from using all concepts
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    
    results_file = ['4ch_phi3Composition_' const_type '_svm_' class_type '.mat'];
    acc_tmp = load([results_dir results_file]);
    
    % Take the values at a specific cost
    cost_param = 1;
    comp_accuracies.(const_type) = permute(mean(acc_tmp.cost_accuracies(cost_param, :, :), 3), [2 3 1]);
    
    % Take the maximum across costs
    comp_accuracies.(const_type) = permute(mean(max(acc_tmp.cost_accuracies, [], 1), 3), [2 3 1]);
    
    cost_acc_all = cat(3, cost_acc_all, permute(acc_tmp.cost_accuracies, [2 3 4 1]));
end

%% Join accuracies together

const_type = const_types{1};

% Single concept results
o_acc_all = o_accuracies.(const_type);

% Full constellation results
comp_acc_all = comp_accuracies.(const_type);

% Join together with big phi accuracies
acc_all = [o_acc_all comp_acc_all p_accuracies];

%% Cost-performance

cost_acc_mean = squeeze(mean(mean(cost_acc_all, 2), 1));
cost_acc_std = squeeze(std(mean(cost_acc_all, 2), [], 1));

plot_measures = [17 16 (1:15)];
measure_offsets = [-0.1 0.1 zeros(1, 15)];

figure;
for measure_c = 1 : length(plot_measures)
    measure = plot_measures(measure_c);
    errorbar((1:size(cost_acc_mean, 2))+measure_offsets(measure_c), cost_acc_mean(measure, :), cost_acc_std(measure, :)); hold on;
end
legend('IIS', 'II');
set(gca, 'XTick', (1:size(cost_acc_mean, 2)), 'XTickLabel', (-20:10:20));
xlabel('cost parameter (2^{x})');

%% Channels vs accuracy
% Each line corresponds to measure-type

channel_sets = nchoosek((1:15), 4);

channel_means = zeros(15, size(acc_all, 2));
channel_stds = zeros(size(channel_means));

for channel = 1 : max(channel_sets(:))
    contains_channel = any(channel_sets == channel, 2);
    channel_means(channel, :) = mean(acc_all(contains_channel, :), 1);
    channel_stds(channel, :) = std(acc_all(contains_channel, :), [], 1);
end

figure;
measure_offsets = [-0.25 -0.15 -0.05 0.05 0.15 0.25];
for measure = 1 : size(channel_means, 2)
    errorbar((1:15)+measure_offsets(measure), channel_means(:, measure), channel_stds(:, measure)); hold on;
end
legend('1', '2', '3', '4', 'all', '\Phi');

%% Channels vs accuracy (2)
% Each line corresponds to a measure-type
% Average across all sets NOT containing the channel

channel_sets = nchoosek((1:15), 4);

channel_means = zeros(15, size(acc_all, 2));
channel_stds = zeros(size(channel_means));

for channel = 1 : max(channel_sets(:))
    contains_channel = ~any(channel_sets == channel, 2);
    channel_means(channel, :) = mean(acc_all(contains_channel, :), 1);
    channel_stds(channel, :) = std(acc_all(contains_channel, :), [], 1);
end

figure;
measure_offsets = [-0.25 -0.15 -0.05 0.05 0.15 0.25];
for measure = 1 : size(channel_means, 2)
    errorbar((1:15)+measure_offsets(measure), channel_means(:, measure), channel_stds(:, measure)); hold on;
end
legend('1', '2', '3', '4', 'all', '\Phi');

%% Channels vs phi
% Each line corresponds to a measure-type
% Average across all sets NOT containing the channel

load('../phi_3/results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM0.mat');
vals = mean(mean(phis{1}.phis(:, :, :, 1) - phis{1}.phis(:, :, :, 2), 2), 3);
%vals = mean(mean(phis{1}.phis(:, :, :, 2), 2), 3);

channel_sets = nchoosek((1:15), 4);

channel_means = zeros(15, size(acc_all, 2));
channel_stds = zeros(size(channel_means));

for channel = 1 : max(channel_sets(:))
    contains_channel = ~any(channel_sets == channel, 2);
    channel_means(channel, :) = mean(vals(contains_channel), 1);
    channel_stds(channel, :) = std(vals(contains_channel), [], 1);
end

figure;
measure_offsets = [-0.25 -0.15 -0.05 0.05 0.15 0.25];
for measure = 1 : size(channel_means, 2)
    errorbar((1:15)+measure_offsets(measure), channel_means(:, measure), channel_stds(:, measure)); hold on;
end
legend('1', '2', '3', '4', 'all', '\Phi');

%% Measure-type vs accuracy
% Each line corresponds to a channel

channel_sets = nchoosek((1:15), 4);

channel_means = zeros(15, size(acc_all, 2));
channel_stds = zeros(size(channel_means));

for channel = 1 : max(channel_sets(:))
    contains_channel = ~any(channel_sets == channel, 2);
    channel_means(channel, :) = mean(acc_all(contains_channel, :), 1);
    channel_stds(channel, :) = std(acc_all(contains_channel, :), [], 1);
end

figure;
channel_offsets = linspace(-0.25, 0.25, max(channel_sets(:)));
for channel = 1 : max(channel_sets(:))
    errorbar((1:size(acc_all, 2))+channel_offsets(channel), channel_means(channel, :), channel_stds(channel, :)); hold on;
end

%% Plot accuracy vs set-center

channel_sets = nchoosek((1:15), 4);

set_means = mean(channel_sets, 2);

%figure; scatter(set_means, acc_all(:, 1));

figure; scatter(set_means, acc_all(:, 6));

%% Plot accuracy vs set-distance

addpath('../');

set_distances = channel_set_distances(channel_sets);

figure; scatter(set_distances, acc_all(:, 6));
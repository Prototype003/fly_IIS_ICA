%% Description

%{

Plots and stats for seeing if there is a correlation between some set
attribute (e.g. mean channel location) and classification accuracy

%}

%% Setup

nChannels = 4;
constellation_type = 'unpart'; % 'part'; 'diff'; 'both'
class_type = 'across';

if strcmp(class_type, 'across')
    var_source_dim = 3; % trial dimension
else % if strcmp(class_type, 'within')
    var_source_dim = 4; % fly dimension
end

%% Load phi values

source_dir = '../phi_3/results/';
source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels' num2str(nChannels) '_globalTPM0.mat'];
tic
disp('loading');
tmp = load([source_dir source_file]);
disp('loaded');
toc

phis = tmp.phis{1};

%% Compute channel set locations (mean of channels)

channel_sets = nchoosek((1:15), 4);

set_means = mean(channel_sets, 2);

%% Preprocess phi values

% Weighted average across states
big_mips = phis.big_mips;
for const_type = 1 : size(phis.big_mips, 2)
    for concept = 1 : size(phis.big_mips, 3)
        big_mips(:, const_type, concept, :, :, :, :) = permute(big_mips(:, const_type, concept, :, :, :, :), [1 4 5 6 7 2 3]) .* single(phis.state_counters);
    end
end
big_mips = permute(sum(big_mips, 1), [2 3 4 5 6 7 1]);
big_mips = big_mips ./ single(sum(phis.state_counters(:, 1, 1, 1, 1))); % State count should be constant across all parameters

% Get desired constellation type
% And specify index-vector (start of each constellation, 0-indexed as any
%   counters will be added to it, i.e. const_starts+counter)
if strcmp(constellation_type, 'diff')
    % Unpartitioned - partitioned
    big_mips = permute(...
        big_mips(1, :, :, :, :, :) - big_mips(2, :, :, :, :, :),...
        [2 3 4 5 6 1]...
        );
    const_starts = [0];
elseif strcmp(constellation_type, 'both')
    % Both unpartitioned and partitioned
    big_mips = permute(...
        cat(2, big_mips(1, :, :, :, :, :), big_mips(2, :, :, :, :, :)),...
        [2 3 4 5 6 1]...
        );
    const_starts = [0 size(phis.big_mips, 3)];
elseif strcmp(constellation_type, 'unpart')
    % Unpartitioned
    big_mips = permute(...
        big_mips(1, :, :, :, :, :),...
        [2 3 4 5 6 1]...
        );
    const_starts = [0];
elseif strcmp(constellation_type, 'part')
    % Unpartitioned
    big_mips = permute(...
        big_mips(2, :, :, :, :, :),...
        [2 3 4 5 6 1]...
        );
    const_starts = [0];
end

%% Correlation (set-mean, big-phi)

aw = zeros(size(phis.phis, 3), 1);
an = zeros(size(aw));
figure;
for fly = 1 : size(phis.phis, 3)
    subplot(3, 5, fly);
    wake_vals = mean(phis.phis(:, :, fly, 1), 2); % mean across trials
    anest_vals = mean(phis.phis(:, :, fly, 2), 2); % mean across trials
    scatter(set_means, wake_vals, 'r.');
    hold on;
    scatter(set_means, anest_vals, 'b.');
    aw(fly) = corr(set_means, wake_vals);
    an(fly) = corr(set_means, anest_vals);
end

% one-sample t-test from 0
[h, p, ci, stats] = ttest(fisher_rz(aw))
[h, p, ci, stats] = ttest(fisher_rz(an))

%% Correlation (set-mean, small-phis)

aw = zeros(size(big_mips, 4), size(big_mips, 1));
an = zeros(size(aw));
for concept = 1 : size(big_mips, 1)
    for fly = 1 : size(big_mips, 4)
        wake_vals = mean(big_mips(concept, :, :, fly, 1), 3); % mean across trials
        anest_vals = mean(big_mips(concept, :, :, fly, 2), 3); % mean across trials
        aw(fly, concept) = corr(set_means, wake_vals(:));
        an(fly, concept) = corr(set_means, anest_vals(:));
    end
end

% one-sample t-tests from 0
clear hs ps
for concept = 1 : size(big_mips, 1)
    [hs(concept, 1), ps(concept, 1)] = ttest(aw(:, concept)); % wake
    [hs(concept, 2), ps(concept, 2)] = ttest(an(:, concept)); % anest
end

%% Load accuracies per concept (and big-phi)

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
    
    acc_repeatMean = c_accuracies.(const_type);
    acc_repeats = permute(max(acc_tmp.cost_accuracies, [], 4), [1 3 2]);
    
end

%% Load accuracies for full IIS

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
    
    acc_repeatMean = cat(2, acc_repeatMean, comp_accuracies.(const_type));
    acc_repeats = cat(2, acc_repeats, permute(max(acc_tmp.cost_accuracies, [], 1), [2 1 3]));
    
end

%% Correlation (set-mean, accuracies) - after averaging across repeats (flies or trials)

rs = zeros(size(acc_repeatMean, 2), 1);
ps = zeros(size(rs));
figure;
for measure = 1 : size(acc_repeatMean, 2)
    subplot(3, 6, measure);
    vals = acc_repeatMean(:, measure);
    scatter(set_means, vals, '.');
    [rs(measure), ps(measure)] = corr(set_means, vals);
end

%% Correlation (set-mean, accuracies) - for each individual repeat (flies or trials)

rs = zeros(size(acc_repeats, 2), 1);
for measure = 1 : size(acc_repeats, 2)
    for repeat = 1 : size(acc_repeats, 3)
        vals = acc_repeats(:, measure, repeat);
        [rs(measure, repeat)] = corr(set_means, vals);
    end
end

% one-sample t-tests from 0
clear hs ps
for measure = 1 : size(acc_repeats, 2)
    [hs(measure), ps(measure)] = ttest(fisher_rz(rs(measure, :)))
end

%%

foo = squeeze(acc_repeats(:, end, :));
rs = zeros(size(foo, 2), 1);
ps = zeros(size(rs));
figure;
for i = 1 : size(foo, 2)
    subplot(3, 5, i);
    scatter(set_means, foo(:, i) + (rand([1365 1])./50), '.');
    [rs(i), ps(i)] = corr(set_means, foo(:, i)+(rand([1365 1])./50));
end

%% Direct comparison of central sets versus peripheral sets versus both

% Central - all channels less than 8
% Peripheral - all channels greater than 8
% Both - have channels in both central and peripheral

middle = median((1:15));

center_set_inds = all(channel_sets < middle, 2);
periph_set_inds = all(channel_sets > middle, 2);
both_set_inds = any(channel_sets < middle, 2) & any(channel_sets > middle, 2);

center_sets = channel_sets(center_set_inds, :);
periph_sets = channel_sets(periph_set_inds, :);
both_sets = channel_sets(both_set_inds, :);

% Get accuracies of each group of channels
center_accs = acc_repeats(center_set_inds, :, :);
periph_accs = acc_repeats(periph_set_inds, :, :);


%% Stats (LME)

% SII values
condition = 2;
values = permute(mean(phis.phis(:, :, :, condition), 2), [1 3 2 4]); % mean across trials

% Accuracy values
measure = 17;
values = permute(acc_repeats(:, measure, :), [1 3 2]); % not averaged across repeats

% Create table
tic;
table_array = zeros(numel(values), 4);
row_counter = 1;
for repeat = 1 : size(values, 2)
    for network = 1 : size(values, 1)
        row = [network repeat mean(channel_sets(network, :)) values(network, repeat)];
        table_array(row_counter, :) = row;
        row_counter = row_counter + 1;
    end
end
toc

acc_table = array2table(table_array, 'VariableNames', {'set', 'repeat', 'set_center', 'x'});
acc_table.set = categorical(acc_table.set);
acc_table.repeat = nominal(acc_table.repeat);

%acc_table.x = zscore(acc_table.x);
%acc_table.set_center = zscore(acc_table.set_center);

%% LME - does include network location explain more variance?

% LME
model_full_spec = 'x ~ set_center + (1|repeat)';
model_full = fitlme(acc_table, model_full_spec);

% Null
model_null_spec = 'x ~ 1 + (1|repeat)';
model_null = fitlme(acc_table, model_null_spec);
compare(model_null, model_full)

%% Set-center LME

feature_type = categorical(4); % 5 = composition; 6 = big-phi
feature_table = acc_table(acc_table.feature_type == feature_type, :);

% LME
model_full_spec = 'accuracy ~ set_center + (1|set)';
model_full = fitlme(feature_table, model_full_spec);

% Null
model_null_spec = 'accuracy ~ 1 + (1|set)';
model_null = fitlme(feature_table, model_null_spec);
compare(model_null, model_full)

[r, p] = corr(feature_table.set_center, feature_table.accuracy)
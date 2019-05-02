%% Figure for Nao grant (2019-04-08)

%{

%}

%% Load fly data

load('../workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat');

nTrials = size(fly_data, 3);
nFlies = size(fly_data, 4);
nConditions = size(fly_data, 5);
tau = 1;

%% Compute pairwise correlations

% Can add pre-processing (e.g. binarisation) here
% Median split
fly_medians = median(fly_data, 1);
fly_medians = repmat(fly_medians, [size(fly_data, 1) 1 1 1 1]);
fly_data = fly_data > fly_medians;

% Correlation between each pair of channels
correlations = zeros(size(fly_data, 2), size(fly_data, 2), size(fly_data, 3), size(fly_data, 4), size(fly_data, 5));
for fly = 1 : size(fly_data, 4)
    for condition = 1 : size(fly_data, 5)
        for trial = 1 : size(fly_data, 3)
            r = corr(fly_data(:, :, trial, fly, condition));
            correlations(:, :, trial, fly, condition) = r;
        end
    end
end

% Convert to HCTSA style matrix
% Order - wake/anest > fly > trial
correlation_features = zeros(nTrials*nFlies*nConditions, nchoosek(size(correlations,1), 2));
upper_triangle = logical(tril(ones(size(correlations,1), size(correlations,2)), -1)); % tril to match order of nchoosek (because indexing is row-wise first
row = 1;
for condition = 1 : size(correlations, 5)
    for fly = 1 : size(correlations, 4)
        for trial = 1 : size(correlations, 3)
            tmp_features = correlations(:, :, trial, fly, condition);
            correlation_features(row, :) = tmp_features(upper_triangle);
            row = row + 1;
        end
    end
end

%% Load phi-3 values

addpath('../figure_code/');

% Load phi-3 results (~30 seconds)
tic;
[phis, phi_string] = phi_load('phi_three', 0, '../');
toc

%% phi-3 to HCSTA style matrix

% Convert to HCTSA style matrix
phi_features = [];
for nChannels = 1 : length(phis)
    feature_mat = zeros(nTrials*nFlies*nConditions, size(phis{nChannels}.phis, 1));
    row = 1;
    for condition = 1 : size(correlations, 5)
        for fly = 1 : size(correlations, 4)
            for trial = 1 : size(correlations, 3)
                feature_mat(row, :) = phis{nChannels}.phis(:, trial, fly, condition, tau)';
                row = row + 1;
            end
        end
    end
    phi_features = cat(2, phi_features, feature_mat);
end

%% Normalise HCTSA style features (across ALL flies and trials)

hctsa_mat = [correlation_features phi_features];
%hctsa_mat = hctsa_mat([(1:8) (1:8)+104], :); % For fly 1

% Normalise each column using min-max method
% (https://stats.stackexchange.com/questions/70801/how-to-normalize-data-to-0-1-range)
feature_mins = repmat(min(hctsa_mat, [], 1), [size(hctsa_mat, 1) 1]);
feature_maxes = repmat(max(hctsa_mat, [], 1), [size(hctsa_mat, 1) 1]);
hctsa_mat = (hctsa_mat - feature_mins) ./ (feature_maxes - feature_mins);

%% Normalise HCTSA style features (within EACH fly)

hctsa_mat = [correlation_features phi_features];

rows = [(1:nTrials) (1:nTrials) + size(hctsa_mat, 1)/2];
for fly = 1 : nFlies
    fly_mat = hctsa_mat(rows, :);
    fly_mins = repmat(min(fly_mat, [], 1), [size(fly_mat, 1) 1]);
    fly_maxes = repmat(max(fly_mat, [], 1), [size(fly_mat, 1) 1]);
    hctsa_mat(rows, :) = (fly_mat - fly_mins) ./ (fly_maxes - fly_mins);
    hctsa_mat(rows, :) = zscore(fly_mat);
    rows = rows + nTrials;
end

%% Normalised feature matrix (figure)

clim = prctile(hctsa_mat(:), [5 95]);

figure; imagesc(hctsa_mat, clim); cbar = colorbar; colormap inferno;
hold on;

% Draw section-lines
section_edges = [105.5 105.5+105 105.5+105+455]; % TODO: make dynamic
for section_edge = section_edges
    line([section_edge section_edge], [0 size(hctsa_mat, 1)+1], 'Color', [0 1 0]);
end

% Draw condition lines
line([0 size(hctsa_mat, 2)+1], [104.5 104.5], 'Color', [0 1 0]);

% Draw fly lines
row = 8.5;
for condition = 1 : nConditions
    for fly = 1 : nFlies
        line([0 size(hctsa_mat, 2)+1], [row row], 'Color', [0 1 0], 'LineStyle', ':'); % Draw line at last trial for the fly
        row = row + nTrials;
    end
end

%set(gca, 'YTick', (1:nTrials:size(hctsa_mat, 1)));
set(gca, 'YTick', [52 52+104], 'YTickLabel', {'wake', 'anest'});

set(gca, 'XTick', [52.5 105+52.5 105+105+(455/2) 105+105+455+(1365/2)],...
    'XTickLabel', {'r_{2ch}', '\Phi^{3.0}_{2ch}', '\Phi^{3.0}_{3ch}', '\Phi^{3.0}_{4ch}'});

ylabel(cbar, 'min-max normed value');
ylabel(cbar, 'zscored value');

%% Classify (nearest-mean leave-one-out) correlation

addpath('../');

% Within-fly classification
correlation_accuracies_w = zeros(nFlies, size(correlation_features, 2));
fly_rows = (1:nTrials);
for fly = 1 : nFlies
    class_data = zeros(nTrials, nConditions, size(correlation_features, 2));
    for condition = 1 : nConditions
        class_data(:, condition, :) = correlation_features(fly_rows+((condition-1)*(size(correlation_features, 1)/nConditions)), :);
    end
    for feature = 1 : size(class_data, 3)
        correlation_accuracies_w(fly, feature) = classify_lol(class_data(:, :, feature));
    end
    fly_rows = fly_rows + nTrials;
end

% Across-fly classification
correlation_accuracies_a = zeros(size(correlation_features, 2), 1);
correlations_trMean = permute(mean(correlations, 3), [1 2 4 5 3]);
channel_pairs = nchoosek((1:15), 2);
for p = 1 : size(channel_pairs, 1)
    pair = channel_pairs(p, :);
    class_data = permute(correlations_trMean(pair(1), pair(2), :, :), [3 4 1 2]);
    correlation_accuracies_a(p) = classify_lol(class_data);
end

%% Load phi-3 classification accuracies

phi_accuracies_w = load('../workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat');
phi_accuracies_w = phi_accuracies_w.accuracies;

phi_accuracies_a = load('../workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification_across1.mat');
phi_accuracies_a = phi_accuracies_a.accuracies;

%% Join correlation accuracies with phi accuracies

% Within
accuracies_w = cell(length(phi_accuracies_w)+1, 1);
accuracies_w{1} = struct(); accuracies_w{1}.accuracies = correlation_accuracies_w';
accuracies_w(2:end) = phi_accuracies_w;

% Across
accuracies_a = cell(length(phi_accuracies_a)+1, 1);
accuracies_a{1} = struct(); accuracies_a{1}.accuracies = correlation_accuracies_a;
accuracies_a(2:end) = phi_accuracies_a;

%% Sort HCTSA-style feature matrix by classification accuracy
% Sort by within-fly classification

feature_start = 1;
for feature_type = 1 : length(accuracies_w)
    [~, order] = sort(mean(accuracies_w{feature_type}.accuracies, 2)); % Assumes values to be sorted to be in 1st dimension
    features = (feature_start:feature_start-1+size(accuracies_w{feature_type}.accuracies, 1));
    feature_mat = hctsa_mat(:, features);
    hctsa_mat(:, features) = feature_mat(:, order);
    feature_start = feature_start + size(feature_mat, 2);
end

%% Sort HCTSA-style feature matrix by classification accuracy
% Sort by across-fly classification

feature_start = 1;
for feature_type = 1 : length(accuracies_a)
    [~, order] = sort(mean(accuracies_a{feature_type}.accuracies, 2)); % Assumes values to be sorted to be in 1st dimension
    features = (feature_start:feature_start-1+size(accuracies_a{feature_type}.accuracies, 1));
    feature_mat = hctsa_mat(:, features);
    hctsa_mat(:, features) = feature_mat(:, order);
    feature_start = feature_start + size(feature_mat, 2);
end

%% Classification histogram (within)
% accuracies by number of features

% Bar histogram
figure;
for feature_type = 1 : length(accuracies_w)
    histogram(mean(accuracies_w{feature_type}.accuracies, 2), 'Normalization', 'pdf'); hold on;
end
legend('r', '2ch', '3ch', '4ch');
xlabel('accuracy'); ylabel('p');

% "Line" histogram
figure;
for feature_type = 1 : length(accuracies_w)
    [N, edges] = histcounts(mean(accuracies_w{feature_type}.accuracies, 2), 'Normalization', 'pdf');
    x=filter([0.5 0.5], 1, edges);
    plot(x(2:end), N); hold on;
end
legend('r', '2ch', '3ch', '4ch');
xlabel('accuracy'); ylabel('p');

% Cumulative distribution
figure;
for feature_type = 1 : length(accuracies_w)
    cdfplot(mean(accuracies_w{feature_type}.accuracies, 2)); hold on;
end
legend('r', '2ch', '3ch', '4ch');
xlabel('accuracy'); ylabel('p');
title('within-fly classification');

%% Classification histogram (across)
% accuracies by number of features

% Bar histogram
figure;
for feature_type = 1 : length(accuracies_a)
    histogram(mean(accuracies_a{feature_type}.accuracies, 2), 'Normalization', 'pdf'); hold on;
end
legend('r', '2ch', '3ch', '4ch');
xlabel('accuracy'); ylabel('p');

% "Line" histogram
figure;
for feature_type = 1 : length(accuracies_a)
    [N, edges] = histcounts(mean(accuracies_a{feature_type}.accuracies, 2), 'Normalization', 'pdf');
    x=filter([0.5 0.5], 1, edges);
    plot(x(2:end), N); hold on;
end
legend('r', '2ch', '3ch', '4ch');
xlabel('accuracy'); ylabel('p');

% Cumulative distribution
figure;
for feature_type = 1 : length(accuracies_a)
    cdfplot(mean(accuracies_a{feature_type}.accuracies, 2)); hold on;
end
legend('r', '2ch', '3ch', '4ch');
xlabel('accuracy'); ylabel('p');
title('across-fly classification');
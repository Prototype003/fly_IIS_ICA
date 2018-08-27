%% Description

%{

Classification across trials for a single channel set, single fly

Need to classify:
    Using 2nd order concepts only
    Using 3rd order concepts only
    Using 4th order concepts only
For
    State dependent (all states x concepts)
    Unweighted mean (sum of concept phis across possible states / number of possible states)
    Weighted mean (weighted mean of concept phis across states)

%}

%% Setup

addpath('C:\Users\this_\Documents\MATLAB\Toolboxes\liblinear-2.20\windows');

addpath('../svm_classification/');

nChannels = 4;
flies = (1:13);
conditions = (1:2);
tau = 4;
trials = (1:8);
global_tpm = 0;
tau_bin = 0;
sample_offset = 0;

% split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f03c2tau4tauOffset0s0939t6.mat
file_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels' num2str(nChannels)];

nStates = 2^nChannels;
nConcepts = 15;

% (state-dependence)x(which concepts to use)
classification_features = cell(1, 5);
for concept_column = 1 : size(classification_features, 2)
    classification_features{1, concept_column}.states = nStates+1; classification_features{1, concept_column}.average_states = 0; % Weighted averaged states
end
for state_row = 1 : size(classification_features, 1)
    classification_features{state_row, 1}.concepts = (1:4); % first order
    classification_features{state_row, 2}.concepts = (5:10); % second order
    classification_features{state_row, 3}.concepts = (11:14); % third order
    classification_features{state_row, 4}.concepts = (15); % fourth order
    classification_features{state_row, 5}.concepts = (5:15); % 2nd-4th orders
end

%% Load

source_dir = 'results/';
source_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM0.mat';
tic
disp('loading');
tmp = load([source_dir source_file]);
disp('loaded');
toc

phis = tmp.phis{1};

%% Preprocess

% Unpartitioned - partitioned
big_mips = permute(...
    phis.big_mips(:, 1, :, :, :, :, :, :) - phis.big_mips(:, 2, :, :, :, :, :, :),...
    [1 3 4 5 6 7 8 2]...
    );
% % Unpartitioned
% big_mips = permute(...
%     phis.big_mips(:, 1, :, :, :, :, :, :),...
%     [1 3 4 5 6 7 8 2]...
%     );

% Weighted average across states
for concept = 1 : size(big_mips, 2)
    %big_mips(:, concept, :, :, :, :) = big_mips(:, concept, :, :, :, :) .* single(phis.state_counters);
    big_mips(:, concept, :, :, :, :) = permute(big_mips(:, concept, :, :, :, :), [1 3 4 5 6 7 2]) .* single(phis.state_counters);
end
big_mips = permute(sum(big_mips, 1), [2 3 4 5 6 7 1]);
big_mips = big_mips ./ single(sum(phis.state_counters(:, 1, 1, 1, 1)));

%% Classify

accuracies = zeros(size(big_mips, 4), size(big_mips, 2), size(classification_features, 2)+1);

cost = 1;

for fly = 1 : size(big_mips, 4) % Roughly 70s per fly
    disp(fly); tic;
    for network = 1 : size(big_mips, 2)
        for concept_features = 1 : size(classification_features, 2)
            
            % Get features (concepts)
            concepts = classification_features{1, concept_features}.concepts;
            features = permute(...
                big_mips(concepts, network, :, fly, :),...
                [3 1 5 2 4]...
                ); % trials x features x conditions (observations x features x classes)
            
            % Classify
            results = svm_lol_liblinear_manual(features, cost);
            accuracies(fly, network, concept_features) = results.accuracy;
            
        end
        
        % Big phi classification
        features = permute(phis.phis(network, :, fly, :), [2 1 4 3]);
        results = svm_lol_liblinear_manual(features, cost);
        accuracies(fly, network, concept_features+1) = results.accuracy;
        
    end
    toc
end

%% Save accuracies

results_dir = 'results/';
results_file = 'composition_classify.mat';
save([results_dir results_file], 'accuracies', 'tau', 'nChannels');

%% Load accuracies

results_dir = 'results/';
results_file = 'composition_classify.mat';
load([results_dir results_file]);

%% Plot

order_labels = {'1st', '2nd', '3rd', '4th', '2/3/4th', '\Phi'};

% Plots average across flies, for all channel sets
figure;
imagesc(squeeze(mean(accuracies, 1)));
c = colorbar; ylabel(c, 'accuracy');
set(gca, 'XTickLabels', order_labels);
xlabel('concept-orders used');
ylabel('channel set');

% Plots average across channel sets, and average + std across flies
figure;
errorbar((1:6), squeeze(mean(mean(accuracies, 2), 1)), squeeze(std(mean(accuracies, 2), [], 1)));
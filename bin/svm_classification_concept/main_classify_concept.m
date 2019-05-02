function [] = main_classify_concept(constellation_type)

%% Description

%{

Classification across trials for all channel sets and flies

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

% 'unpart': unpartitioned; 'part': partitioned
% 'diff': unpartitioned-partitioned; 'both': unpartitioned AND partitioned
%constellation_type = 'unpart';

addpath('C:\Users\this_\Documents\MATLAB\Toolboxes\liblinear-2.20\windows');

addpath('../svm_classification/');

results_location = 'results/';

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
if strcmp(constellation_type, 'both') % 2x sets of concept values
    for state_row = 1 : size(classification_features, 1)
        classification_features{state_row, 1}.concepts = [(1:4) (1:4)+15]; % first order
        classification_features{state_row, 2}.concepts = [(5:10) (5:10)+15]; % second order
        classification_features{state_row, 3}.concepts = [(11:14) (11:14)+15]; % third order
        classification_features{state_row, 4}.concepts = [(15) (15)+15]; % fourth order
        %classification_features{state_row, 5}.concepts = [(5:15) (5:15)+15]; % 2nd-4th orders
        %classification_features{state_row, 5}.concepts = [(1:15) (1:15)+15]; % 1st-4th orders
    end
else % only one set of concept values
    for state_row = 1 : size(classification_features, 1)
        classification_features{state_row, 1}.concepts = (1:4); % first order
        classification_features{state_row, 2}.concepts = (5:10); % second order
        classification_features{state_row, 3}.concepts = (11:14); % third order
        classification_features{state_row, 4}.concepts = (15); % fourth order
        %classification_features{state_row, 5}.concepts = (5:15); % 2nd-4th orders
        %classification_features{state_row, 5}.concepts = (1:15); % 1st-4th orders
    end
end

%% Load

source_dir = '../phi_3/results/';
source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels' num2str(nChannels) '_globalTPM0.mat'];
tic
disp('loading');
tmp = load([source_dir source_file]);
disp('loaded');
toc

phis = tmp.phis{1};

%% Preprocess

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

%% Classify across flies (parallel)
% ~12 seconds per network with cost search (-20:10:20)

class_type = 'across';
results_file = [num2str(nChannels) 'ch_phi3Concept_' constellation_type '_svm_' class_type '.mat'];

cost_powers = (-20:10:20);
costs = 2 .^ cost_powers;

% networks x concepts+Phi x costs
cost_accuracies = zeros(size(big_mips, 2), nConcepts+1, length(costs));

% Broadcast variables
phis_p = parallel.pool.Constant(phis.phis);
big_mips_p = parallel.pool.Constant(big_mips);
const_starts_p = parallel.pool.Constant(const_starts);
costs_p = parallel.pool.Constant(costs);

parfor network = 1 : size(big_mips, 2)
    disp(network); tic;
    
    nConcepts = size(big_mips_p.Value, 1) / length(const_starts_p.Value);
    
    for cost_counter = 1 : length(costs_p.Value)
        cost = costs_p.Value(cost_counter);
        
        accuracies = zeros(nConcepts+1, length(costs_p.Value));
        
        for concept = 1 : nConcepts
            concept_inds = const_starts_p.Value + concept;
            
            % Get features (concepts), average across trials
            features = permute(...
                mean(big_mips_p.Value(concept_inds, network, :, :, :), 3),...
                [4 1 5 2 3]...
                ); % flies x features x conditions (observations x features x classes)
            
            % Classify
            results = svm_lol_liblinear_manual(features, cost);
            accuracies(concept, cost_counter) = results.accuracy;
            
        end
        
        % Big phi classification
        features = permute(mean(phis_p.Value(network, :, :, :), 2), [3 1 4 2]);
        results = svm_lol_liblinear_manual(features, cost);
        accuracies(nConcepts+1, cost_counter) = results.accuracy;
        cost_accuracies(network, :, :) = accuracies;
        
    end
   
    toc
end

%% Save accuracies

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs', 'nChannels', 'tau');

disp('saved across');

end
%% Description

%{
Finds correlation between TPM divergence and classification accuracy
%}

%% Setup

const_types = {'unpart', 'part', 'both'};
const_types = const_types(1);
class_type = 'within';

%%

load('results/tpm_divergence.mat');

%% Load accuracies

c_accuracies = struct(); % accuracies for each individual concept

results_dir = '../svm_classification_concept/results/';
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    results_file = ['4ch_phi3Concept_' const_type '_svm_' class_type '.mat'];
    acc_tmp = load([results_dir results_file]);
    
    % Take the values at a specific cost
    cost_param = 1;
    c_accuracies.(const_type) = acc_tmp.cost_accuracies(:, :, :, cost_param);
    
    % Take the maximum across costs
    %c_accuracies.(const_type) = max(acc_tmp.cost_accuracies, [], 3);
    
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
        %o_accuracies.(const_type)(:, order_c) = max(c_accuracies.(const_type)(:, concept_indices), [], 2);
        
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
    comp_accuracies.(const_type) = acc_tmp.cost_accuracies(cost_param, :, :);
    
    % Take the maximum across costs
    %comp_accuracies.(const_type) = max(acc_tmp.cost_accuracies, [], 1);
end

%% Correlation

divs = permute(mean(tpm_divergence, 2), [1 3 4 2]); % nets x flies x conds

cond = 1;

figure;
a = divs(:, :, cond);
b = c_accuracies.(const_type)(:, :, 15);
scatter(a(:), b(:), '.');
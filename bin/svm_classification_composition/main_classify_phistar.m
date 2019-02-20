%% Description

%{

Classification across trials for all channel sets and flies

Using big phi from 2- and 3-channels as compositional structure

%}

%% SETUP

gaussian = 1;
if gaussian == 1
    load_string = 'phi_star_gaussian';
    file_string = 'phiStarGaussianComposition';
else
    load_string = 'phi_star';
    file_string = 'phiStarComposition';
end

nChannels = 4;
tau = 1;
global_tpm = 0;

% bin directory location
bin_location = '../';

addpath('../figure_code/'); % phi-three loading function is here
addpath('C:\Users\this_\Documents\MATLAB\Toolboxes\liblinear-2.20\windows');
addpath('../svm_classification/');

results_location = 'results/';

%% Load phi values

[phi_values, measure_strings{1}] = phi_load(load_string, global_tpm, bin_location);

values_all = phi_values{nChannels - 1}.phis(:, :, :, :, tau);
channel_sets = phi_values{nChannels - 1}.channel_sets;

%% Get "compositions" for each channel set

% Composition consists of big-phi of subsets

% Number of concepts
nConcepts = 0;
for c_order = 2 : nChannels % start at two because there are no phi values for 1-channel
    nConcepts = nConcepts + nchoosek(nChannels, c_order);
end

big_mips = zeros(nConcepts, size(channel_sets, 1), size(values_all, 2), size(values_all, 3), size(values_all, 4)); % comp-phis x sets x trials x flies x conditions

for network = 1 : size(channel_sets, 1)
    
    concept_counter = 1;
    
    for c_order = 2 : nChannels 
        
        % Find sets which are subsets
        subset_inds = all(ismember(phi_values{c_order-1}.channel_sets, channel_sets(network, :)), 2);
        
        % Get corresponding phi values for subsets
        big_mips((concept_counter:concept_counter+sum(subset_inds)-1), network, :, :, :) = permute(...
            phi_values{c_order-1}.phis(subset_inds, :, :, :, tau),...
            [1 5 2 3 4]);
        
        concept_counter = concept_counter + sum(subset_inds);
    end
end

%% Classify within fly
% ~10 hours for all flies, with cost search (when using just 1 cpu)
% ~300 seconds for 1 fly, with cost search (when using 4 cpus)

class_type = 'within';
results_file = [file_string '_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(big_mips, 2), size(big_mips, 4));

big_mips_par = parallel.pool.Constant(big_mips);
costs_par = parallel.pool.Constant(costs);

for fly = 1 : size(big_mips, 4) % Roughly 70s per fly, cost-level
    disp(['fly ' num2str(fly)]); tic;
    
    parfor network = 1 : size(big_mips, 2)
        disp(network);
        
        features = permute(big_mips_par.Value(:, network, :, fly, :), [3 1 5 2 4]); % trials x comp-phis x conditions
        accuracies = zeros(size(costs_par.Value));
        
        for cost_counter = 1 : length(costs_par.Value)
            cost = costs_par.Value(cost_counter);
            results = svm_lol_liblinear_manual(features, cost); % observations x features x classes
            accuracies(cost_counter) = results.accuracy;
        end
        
        cost_accuracies(:, network, fly) = accuracies;
        
    end
    
    toc
    
end

accuracy = zeros(size(cost_accuracies, 2), size(cost_accuracies, 3));
for fly = 1 : size(cost_accuracies, 3)
    for network = 1 : size(cost_accuracies, 2)
        accuracy(network, fly) = max(cost_accuracies(:, network, fly));
    end
end

%% Save accuracies

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs', 'nChannels');

disp('saved within');

%% Classify across flies
% ~40 minutes

class_type = 'across';
results_file = [file_string '_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(big_mips, 2));

big_mips_par = parallel.pool.Constant(big_mips);
costs_par = parallel.pool.Constant(costs);

tic;
parfor network = 1 : size(big_mips, 2)
    disp(network); tic;
    
    features = permute(mean(big_mips_par.Value(:, network, :, :, :), 3), [4 1 5 2 3]); % flies x comp-phis x conditions
    accuracies = zeros(size(costs_par.Value));
    
    for cost_counter = 1 : length(costs_par.Value)
        cost = costs_par.Value(cost_counter);
        
        results = svm_lol_liblinear_manual(features, cost);
        accuracies(cost_counter) = results.accuracy;
        
    end
    
    cost_accuracies(:, network) = accuracies;

end
toc
accuracy = max(cost_accuracies, [], 1);

%% Save accuracies

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs', 'nChannels');

disp('saved across');

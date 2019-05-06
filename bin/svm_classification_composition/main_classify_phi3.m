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
constellation_type = 'both';

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

if strcmp(constellation_type, 'diff')
    % Unpartitioned - partitioned
    big_mips = permute(...
        phis.big_mips(:, 1, :, :, :, :, :, :) - phis.big_mips(:, 2, :, :, :, :, :, :),...
        [1 3 4 5 6 7 8 2]...
        );
elseif strcmp(constellation_type, 'both')
    % Both unpartitioned and partitioned
    big_mips = permute(...
        cat(3, phis.big_mips(:, 1, :, :, :, :, :, :), phis.big_mips(:, 2, :, :, :, :, :, :)),...
        [1 3 4 5 6 7 8 2]...
        );
elseif strcmp(constellation_type, 'unpart')
    % Unpartitioned
    big_mips = permute(...
        phis.big_mips(:, 1, :, :, :, :, :, :),...
        [1 3 4 5 6 7 8 2]...
        );
elseif strcmp(constellation_type, 'part')
    % Unpartitioned
    big_mips = permute(...
        phis.big_mips(:, 2, :, :, :, :, :, :),...
        [1 3 4 5 6 7 8 2]...
        );
end

% Weighted average across states
for concept = 1 : size(big_mips, 2)
    %big_mips(:, concept, :, :, :, :) = big_mips(:, concept, :, :, :, :) .* single(phis.state_counters);
    big_mips(:, concept, :, :, :, :) = permute(big_mips(:, concept, :, :, :, :), [1 3 4 5 6 7 2]) .* single(phis.state_counters);
end
big_mips = permute(sum(big_mips, 1), [2 3 4 5 6 7 1]);
big_mips = big_mips ./ single(sum(phis.state_counters(:, 1, 1, 1, 1)));

%% Classify within fly
% ~10 hours for all flies, with cost search (when using just 1 cpu)
% ~300 seconds for 1 fly, with cost search (when using 4 cpus)

% class_type = 'within';
% results_file = [num2str(nChannels) 'ch_phi3Composition_' constellation_type '_svm_' class_type '.mat'];
% 
% cost_powers = (-20:20);%0;%(-20:20);
% costs = 2 .^ cost_powers;
% cost_accuracies = zeros(length(costs), size(big_mips, 2), size(big_mips, 4));
% 
% big_mips_par = parallel.pool.Constant(big_mips);
% costs_par = parallel.pool.Constant(costs);
% 
% for fly = 1 : size(big_mips, 4) % Roughly 70s per fly, cost-level
%     disp(['fly ' num2str(fly)]); tic;
%     
%     parfor network = 1 : size(big_mips, 2)
%         disp(network);
%         
%         features = permute(big_mips_par.Value(:, network, :, fly, :), [3 1 5 2 4]); % trials x comp-phis x conditions
%         accuracies = zeros(size(costs_par.Value));
%         
%         for cost_counter = 1 : length(costs_par.Value)
%             cost = costs_par.Value(cost_counter);
%             results = svm_lol_liblinear_manual(features, cost); % observations x features x classes
%             accuracies(cost_counter) = results.accuracy;
%         end
%         
%         cost_accuracies(:, network, fly) = accuracies;
%         
%     end
%     
%     toc
%     
% end
% 
% accuracy = zeros(size(cost_accuracies, 2), size(cost_accuracies, 3));
% for fly = 1 : size(cost_accuracies, 3)
%     for network = 1 : size(cost_accuracies, 2)
%         accuracy(network, fly) = max(cost_accuracies(:, network, fly));
%     end
% end
% 
% %% Save accuracies
% 
% save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs', 'nChannels', 'tau');
% 
% disp('saved within');

%% Classify across flies
% ~40 minutes

class_type = 'across';
results_file = [num2str(nChannels) 'ch_phi3Composition_' constellation_type '_svm_' class_type '.mat'];

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

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs', 'nChannels', 'tau');

disp('saved across');

%% Load accuracies

% results_dir = 'results/';
% results_file = ['phi3Composition_' constellation_type '_svm_' class_type '.mat'];
% load([results_dir results_file]);

%% Plot

% flies = (1:13);
% 
% order_labels = {'1st', '2nd', '3rd', '4th', '1/2/3/4th', '\Phi'};
% 
% % Plots average across flies, for all channel sets
% figure;
% imagesc(squeeze(mean(accuracies(flies, :, 1:6), 1)));
% c = colorbar; ylabel(c, 'accuracy');
% set(gca, 'XTickLabels', order_labels);
% xlabel('concept-orders used');
% ylabel('channel set');
% 
% % Plots average across channel sets, and average + std across flies
% figure;
% errorbar((1:6), squeeze(mean(mean(accuracies(flies, :, :), 2), 1)), squeeze(std(mean(accuracies(flies, :, :), 2), [], 1) / size(accuracies, 1)));
% %errorbar((1:6), squeeze(mean(mean(accuracies, 2), 1)), squeeze(std(mean(accuracies, 2), [], 1)));
% set(gca, 'XTickLabels', order_labels);
% xlabel('concept-orders used');
% ylabel('accuracy');
% xlim([0.5 6.5]);
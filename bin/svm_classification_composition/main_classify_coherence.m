%% DESCRIPTION

%{

Classifies awake/anest using 105 features - coherence of 15 channels (at lowest frequency)

Across flies classification
Within flies classification

%}

%% SETUP

nChannels = 4;
networks = nchoosek((1:15), 4); % 15 channels in total
pairs = nchoosek((1:15), 2); % coherence is always pairwise

freq_range = (1:42); % (1:42) = 0-5Hz; Vector of which frequency bins to use, we want to use the lowest frequency bin
freq_range = (43:165); % 10-20Hz

% bin directory location
bin_dir = '../';

results_location = 'results/';

addpath('../svm_classification/');

%% Load coherence

coherencies = load_coherence(bin_dir);

% Mean across frequency bins
coherencies = mean(coherencies(freq_range, :, :, :, :), 1);
coherencies = permute(coherencies, [2 3 4 5 1]); % trials x channels x conditions x flies

%% Determine which channel pairs correspond to which channel sets a priori

% nchannels choose 2 (coherence is always pairwise)
network_pairs = zeros(size(networks, 1), nchoosek(size(networks, 2), 2));

for network_counter = 1 : size(networks, 1)
    network = networks(network_counter, :);
    
    % Every possible pair for this set of channels
    matching_pairs = nchoosek(network, 2);
    
    % Find indices of each pair
    for pair_counter = 1 : size(matching_pairs, 1)
        pair = matching_pairs(pair_counter, :);
        [match, index] = ismember(pair, pairs, 'rows'); % there should only be one match
        network_pairs(network_counter, pair_counter) = index;
    end
    
end

%% Classify across flies

class_type = 'across';
results_file = ['coherence_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(size(costs));

% These create a parallel pool
coherencies_par = parallel.pool.Constant(coherencies);
costs_par = parallel.pool.Constant(costs);
poolobj = gcp;
addAttachedFiles(poolobj,{'../svm_classification/train.mexw64','../svm_classification/train.mexw64'})

tic;
parfor network_counter = 1 : size(networks, 1)
    %tic;
    network = network_pairs(network_counter, :);
    disp(network_counter);
    
    % Average across trials
    values = permute(mean(coherencies_par.Value(:, network, :, :), 1), [4 2 3 1]); % flies x pairs x conditions
    
    accuracies = zeros(size(costs_par.Value));
    
    for cost_counter = 1 : length(costs_par.Value)
        cost = costs_par.Value(cost_counter);
        results = svm_lol_liblinear_manual(values, cost);
        accuracies(cost_counter) = results.accuracy;
    end
    
    cost_accuracies(:, network_counter) = accuracies;
    
    %toc
end

toc

accuracy = max(cost_accuracies, [], 1);

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved across');

%% Classify within flies

class_type = 'within';
results_file = ['coherence_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(coherencies, 4));

% These create a parallel pool
coherencies_par = parallel.pool.Constant(coherencies);
costs_par = parallel.pool.Constant(costs);
poolobj = gcp;
addAttachedFiles(poolobj,{'../svm_classification/train.mexw64','../svm_classification/predict.mexw64'})

for fly = 1 : size(coherencies, 4)
    disp(['fly ' num2str(fly)]);
    
    tic;
    
    parfor network_counter = 1 : size(networks, 1)
        disp(network_counter);
        network = network_pairs(network_counter, :);
        
        values = coherencies_par.Value(:, network, :, fly); % trials x channels x conditions
        accuracies = zeros(size(costs_par.Value));
        
        for cost_counter = 1 : length(costs_par.Value)
            cost = costs_par.Value(cost_counter);
            results = svm_lol_liblinear_manual(values, cost);
            accuracies(cost_counter) = results.accuracy;
        end
        
        cost_accuracies(:, network_counter, fly) = accuracies;
        
    end
    
    toc
    
end

% accuracy = zeros(size(cost_accuracies, 2), size(cost_accuracies, 3));
% for fly = 1 : size(powers, 4)
%     for network = 1 : size(cost_accuracies, 2)
%         accuracy(network, fly) = max(cost_accuracies(:, network, fly));
%     end
% end

accuracy = max(cost_accuracies, [], 1);

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved within');
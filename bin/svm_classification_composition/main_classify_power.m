%% DESCRIPTION

%{

Classifies awake/anest using 15 features - power of 15 channels (at lowest frequency)

Across flies classification
Within flies classification

%}

%% SETUP

nChannels = 4;
networks = nchoosek((1:15), 4); % 15 channels in total

freq_range = (1:42); % (1:42) = 0-5Hz; Vector of which frequency bins to use, we want to use the lowest frequency bin
freq_range = (43:165); % 10-20Hz

% bin directory location
bin_dir = '../';

results_location = 'results/';

addpath('../svm_classification/');

%% Load power

powers = load_power(bin_dir);

% Mean across frequency bins
powers = mean(powers(freq_range, :, :, :, :), 1);
powers = permute(powers, [2 3 4 5 1]); % trials x channels x conditions x flies

%% Classify across flies

class_type = 'across';
results_file = ['power_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(networks, 1));

for network_counter = 1 : size(networks, 1)
    tic;
    network = networks(network_counter, :);
    
    % Average across trials
    values = permute(mean(powers(:, network, :, :), 1), [4 2 3 1]); % flies x channels x conditions
    
    for cost_counter = 1 : length(costs)
        cost = costs(cost_counter);
        results = svm_lol_liblinear_manual(values, cost);
        cost_accuracies(cost_counter, network_counter) = results.accuracy;
    end
    toc
end

accuracy = max(cost_accuracies, [], 1);

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved across');

%% Classify within flies

class_type = 'within';
results_file = ['power_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(networks, 1), size(powers, 4));

for fly = 1 : size(powers, 4)
    disp(num2str(fly));
    
    for network_counter = 1 : size(networks, 1)
        
        values = powers(:, :, :, fly); % trials x channels x conditions
        
        for cost_counter = 1 : length(costs)
            cost = costs(cost_counter);
            results = svm_lol_liblinear_manual(values, cost);
            cost_accuracies(cost_counter, network_counter, fly) = results.accuracy;
        end
        
    end
    
end

accuracy = zeros(size(cost_accuracies, 2), size(cost_accuracies, 3));
for fly = 1 : size(powers, 4)
    for network = 1 : size(cost_accuracies, 2)
        accuracy(network, fly) = max(cost_accuracies(:, network, fly));
    end
end

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved within');
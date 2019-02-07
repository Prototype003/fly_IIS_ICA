%% DESCRIPTION

%{

Classifies awake/anest using 105 features - coherence of 15 channels (at lowest frequency)

Across flies classification
Within flies classification

%}

%% SETUP

class_type = 'within'; % 'across' or 'within'

freq_range = (1:42); % (1:42) = 0-5Hz; Vector of which frequency bins to use, we want to use the lowest frequency bin
freq_range = (43:165); % 10-20Hz

% bin directory location
bin_dir = '../';

results_location = 'results/';

addpath('../svm_classification/');

%% Load power

coherencies = load_coherence(bin_dir);

% Mean across frequency bins
coherencies = mean(coherencies(freq_range, :, :, :, :), 1);
coherencies = permute(coherencies, [2 3 4 5 1]); % trials x channels x conditions x flies

%% Classify across flies

class_type = 'across';
results_file = ['coherence_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(size(costs));

% Average across trials
values = permute(mean(coherencies, 1), [4 2 3 1]); % flies x channels x conditions

for cost_counter = 1 : length(costs)
    cost = costs(cost_counter);
    results = svm_lol_liblinear_manual(values, cost);
    cost_accuracies(cost_counter) = results.accuracy;
end

accuracy = max(cost_accuracies)

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved across');

%% Classify within flies

class_type = 'within';
results_file = ['coherence_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(coherencies, 4));

for fly = 1 : size(coherencies, 4)
    disp(num2str(fly));
    values = coherencies(:, :, :, fly); % trials x channels x conditions
    
    for cost_counter = 1 : length(costs) % ~60 seconds
        cost = costs(cost_counter);
        results = svm_lol_liblinear_manual(values, cost);
        cost_accuracies(cost_counter, fly) = results.accuracy;
    end
    
end

accuracy = zeros(size(cost_accuracies, 2), 1);
for fly = 1 : size(cost_accuracies, 2)
    accuracy(fly) = max(cost_accuracies(:, fly));
end

accuracy_mean = mean(accuracy)

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved within');
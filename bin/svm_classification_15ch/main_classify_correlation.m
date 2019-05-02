%% DESCRIPTION

%{

Classifies awake/anest using 105 features - coherence of 15 channels (at lowest frequency)

Across flies classification
Within flies classification

%}

%% SETUP

class_type = 'within'; % 'across' or 'within'

% bin directory location
bin_dir = '../';

results_location = 'results/';

addpath('../svm_classification/');

%% Compute correlations
% Matrix should be trials x pairs x conditions x flies

load('../workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat');

nTrials = size(fly_data, 3);
nFlies = size(fly_data, 4);
nConditions = size(fly_data, 5);
pairs = nchoosek((1:size(fly_data, 2)), 2);

% Can add pre-processing (e.g. binarisation) here
% Median split
fly_medians = median(fly_data, 1);
fly_medians = repmat(fly_medians, [size(fly_data, 1) 1 1 1 1]);
fly_data = fly_data > fly_medians;

% Correlation between each pair of channels
upper_triangle = logical(tril(ones(size(fly_data,2), size(fly_data,2)), -1)); % tril to match order of nchoosek (because indexing is row-wise first
correlations = zeros(nTrials, size(pairs, 1), nConditions, nFlies);
for fly = 1 : size(fly_data, 4)
    for condition = 1 : size(fly_data, 5)
        for trial = 1 : size(fly_data, 3)
            r = corr(fly_data(:, :, trial, fly, condition));
            correlations(trial, :, condition, fly) = r(upper_triangle);
        end
    end
end

%% Classify across flies

class_type = 'across';
results_file = ['medianSplit_correlation_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(size(costs));

% Average across trials
values = permute(mean(correlations, 1), [4 2 3 1]); % flies x channels x conditions

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
results_file = ['medianSplit_correlation_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(correlations, 4));

for fly = 1 : size(correlations, 4)
    disp(num2str(fly));
    values = correlations(:, :, :, fly); % trials x channels x conditions
    
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
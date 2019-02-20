%% DESCRIPTION

%{

Classifies awake/anest using 15 features - 15 phis for 15 sets of 4 channels (spanning 15 channels)

Across flies classification
Within flies classification

%}

%% SETUP

gaussian = 0;
if gaussian == 1
    load_string = 'phi_star_gaussian';
    file_string = 'phiStarGaussian';
else
    load_string = 'phi_star';
    file_string = 'phiStar';
end

tau = 1;
global_tpm = 0;

% bin directory location
bin_location = '../';

addpath('../figure_code/'); % phi-three loading function is here
addpath('C:\Users\this_\Documents\MATLAB\Toolboxes\liblinear-2.20\windows');
addpath('../svm_classification/');

results_location = 'results/';

%% Load phi values

[phi_threes, measure_strings{1}] = phi_load(load_string, global_tpm, bin_location);

values_all = phi_threes{3}.phis(:, :, :, :, tau);
channel_sets = phi_threes{3}.channel_sets;

%% Classify across flies

class_type = 'across';
results_file = [file_string '_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(size(costs));

% Average across trials
values = permute(mean(values_all, 2), [3 1 4 2]); % flies x features x conditions

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
results_file = [file_string '_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(values_all, 3));

for fly = 1 : size(values_all, 3)
    disp(num2str(fly)); tic;
    values = permute(values_all(:, :, fly, :), [2 1 4 3]); % trials x features x conditions
    
    for cost_counter = 1 : length(costs) % ~60 seconds
        cost = costs(cost_counter);
        results = svm_lol_liblinear_manual(values, cost);
        cost_accuracies(cost_counter, fly) = results.accuracy;
    end
    toc
end

accuracy = zeros(size(cost_accuracies, 2), 1);
for fly = 1 : size(cost_accuracies, 2)
    accuracy(fly) = max(cost_accuracies(:, fly));
end

accuracy_mean = mean(accuracy)

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved within');

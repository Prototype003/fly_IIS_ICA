%% Description

%{

Treat 2- and 3-channel phi as small-phi values

%}

%% SETUP

tau = 1;
global_tpm = 0;

% bin directory location
bin_location = '../';

addpath('../figure_code/'); % phi-three loading function is here
addpath('C:\Users\this_\Documents\MATLAB\Toolboxes\liblinear-2.20\windows');
addpath('../svm_classification/');

results_location = 'results/';

%% Load phi values

[phi_values, measure_strings{1}] = phi_load('phi_star_gaussian', global_tpm, bin_location);

values_all = phi_values{3}.phis(:, :, :, :, tau);
channel_sets = phi_values{3}.channel_sets;

%% Concatenate ALL channel sets

compositions = [];
for nChannels = 1 : length(phi_values)
    
    compositions = cat(1, compositions, phi_values{nChannels}.phis);
    
end

%% Classify across flies

class_type = 'across';
results_file = ['phiStarGaussianComposition_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(size(costs));

% Average across trials
values = permute(mean(compositions(:, :, :, :, tau), 2), [3 1 4 2]); % flies x features x conditions

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
results_file = ['phiStarGaussianComposition_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(compositions, 3));

for fly = 1 : size(compositions, 3)
    disp(num2str(fly)); tic;
    values = permute(compositions(:, :, fly, :, tau), [2 1 4 3]); % trials x features x conditions
    
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

%% MDS
comp_dims = size(compositions(:, :, :, :, tau));

rep_mat = zeros(comp_dims(1), prod(comp_dims(2:end)));

col = 1;
for condition = 1 : comp_dims(4)
    for fly = 1 : comp_dims(3)
        for trial = 1 : comp_dims(2)
            rep_mat(:, col) = compositions(:, trial, fly, condition, tau);
            col = col+1;
        end
    end
end

rep_dis_mat = 1 - corr(rep_mat);

coords = cmdscale(rep_dis_mat);

figure;
scatter(coords(1:104, 1), coords(1:104, 2)); hold on;
scatter(coords((105:end), 1), coords((105:end), 2));
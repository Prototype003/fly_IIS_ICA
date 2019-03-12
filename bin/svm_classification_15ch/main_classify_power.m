%% DESCRIPTION

%{

Classifies awake/anest using 15 features - power of 15 channels (at lowest frequency)

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

powers = load_power(bin_dir);

% Mean across frequency bins
powers = mean(powers(freq_range, :, :, :, :), 1);
powers = permute(powers, [2 3 4 5 1]); % trials x channels x conditions x flies

%% Classify across flies

class_type = 'across';
results_file = ['medianSplit_power_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(size(costs));

% Average across trials
values = permute(mean(powers, 1), [4 2 3 1]); % flies x channels x conditions

for cost_counter = 1 : length(costs)
    cost = costs(cost_counter);
    results = svm_lol_liblinear_manual(values, cost);
    cost_accuracies(cost_counter) = results.accuracy;
end

accuracy = max(cost_accuracies);

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved across');

%% Classify within flies

class_type = 'within';
results_file = ['medianSplit_power_svm_' class_type '.mat'];

cost_powers = (-20:20);%0;%(-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(length(costs), size(powers, 4));

for fly = 1 : size(powers, 4)
    disp(num2str(fly));
    values = powers(:, :, :, fly); % trials x channels x conditions
    
    for cost_counter = 1 : length(costs) % ~60 seconds
        cost = costs(cost_counter);
        results = svm_lol_liblinear_manual(values, cost);
        cost_accuracies(cost_counter, fly) = results.accuracy;
    end
    
end

accuracy = zeros(size(cost_accuracies, 2), 1);
for fly = 1 : size(powers, 4)
    accuracy(fly) = max(cost_accuracies(:, fly));
end

accuracy_mean = mean(accuracy)

%% Save

save([results_location results_file], 'accuracy', 'cost_accuracies', 'costs');

disp('saved within');

%% MDS using each power value

comp_vals = permute(powers, [2 1 4 3]);
comp_dims = size(comp_vals);

rep_mat = zeros(comp_dims(1), prod(comp_dims(2:end)));

col = 1;
for condition = 1 : comp_dims(4)
    for fly = 1 : comp_dims(3)
        for trial = 1 : comp_dims(2)
            rep_mat(:, col) = comp_vals(:, trial, fly, condition);
            col = col+1;
        end
    end
end

rep_dis_mat = 1 - corr(rep_mat);

coords = cmdscale(rep_dis_mat);

figure;
fly_colours = repmat((1:comp_dims(3)), comp_dims(2));
fly_colours = fly_colours(:); fly_colours = fly_colours(1:prod(comp_dims(2:3)));
scatter(coords(1:prod(comp_dims(2:3)), 1), coords(1:prod(comp_dims(2:3)), 2), [], fly_colours, 'o'); hold on;
scatter(coords((prod(comp_dims(2:3))+1:end), 1), coords((prod(comp_dims(2:3))+1:end), 2), [], fly_colours, 'x');


figure; plot(coords((1:prod(comp_dims(2:3))), 1)-coords((prod(comp_dims(2:3))+1:end), 1));
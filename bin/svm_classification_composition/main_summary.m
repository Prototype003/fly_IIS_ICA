%% Description

%{

Summarise SVM classification results

%}

%% Setup

results_location = 'results/';
measures = {'power', 'coherence', 'phi3Composition_unpart'};
measure_labels = {'P', '\phicomp'};

clim = [0 1];

%% Within-fly classification

results = cell(size(measures));
results_mat = []; % set x fly x measure

class_type = 'within';

for measure = 1 : length(measures)
    filename = [measures{measure} '_svm_' class_type '.mat'];
    results{measure} = load([results_location filename]);
    results{measure}.accuracy = squeeze(results{measure}.accuracy); % Power/coherence have a leading singleton dimension
    results_mat = cat(3, results_mat, results{measure}.accuracy);
end

% Average accuracies across flies
nFlies = size(results_mat, 2);
accuracy = permute(mean(results_mat, 2), [1 3 2]);

% Plot accuracy of each channel set, for each measure
figure; imagesc(accuracy, clim); colorbar;

%% Across-fly classification

results = cell(size(measures));
results_mat = []; % set x measure

class_type = 'across';

for measure = 1 : length(measures)
    filename = [measures{measure} '_svm_' class_type '.mat'];
    results{measure} = load([results_location filename]);
    results{measure}.accuracy = permute(results{measure}.accuracy, [2 1]);
    results_mat = cat(2, results_mat, results{measure}.accuracy);
end

accuracy = results_mat;

% Plot accuracy of each channel set, for each measure
figure; imagesc(accuracy, clim); colorbar;
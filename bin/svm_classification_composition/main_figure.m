%% Description

%{

Make summary figure

Average across sets or across flies?

Correlation between measures - does phi perform better when power also
performs better?

%}

%% Setup

results_location = 'results/';
measures = {'power', 'coherence', 'phi3Composition_unpart', 'phiStarComposition', 'phiStarGaussianComposition'};
measure_labels = {'P', 'C', '\phi3comp', '\Phi*comp', '\Phi*gcomp'};

clim = [0 1];

%% Load within-fly results

results = cell(size(measures));
results_mat = []; % set x fly x measure

class_type = 'within';

for measure = 1 : length(measures)
    filename = [measures{measure} '_svm_' class_type '.mat'];
    results{measure} = load([results_location filename]);
    results{measure}.accuracy = squeeze(results{measure}.accuracy); % Power/coherence have a leading singleton dimension
    results_mat = cat(3, results_mat, results{measure}.accuracy);
end

% Average across flies
results_all = permute(mean(results_mat, 2), [1 3 2]);

%% Load across-fly results

results = cell(size(measures));
results_mat = []; % set x measure

class_type = 'across';

for measure = 1 : length(measures)
    filename = [measures{measure} '_svm_' class_type '.mat'];
    results{measure} = load([results_location filename]);
    results{measure}.accuracy = permute(results{measure}.accuracy, [2 1]);
    results_mat = cat(2, results_mat, results{measure}.accuracy);
end

% Join with other results
results_all = cat(3, results_all, results_mat);

%% Make figure


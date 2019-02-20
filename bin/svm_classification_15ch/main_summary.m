%% Description

%{

Summarise power and phi-3 SVM classification results

%}

%% Setup

results_location = 'results/';
measures = {'power', 'coherence', 'phi3', 'phi3CompositionOverlapping', 'PhiStar', 'PhiStarComposition', 'PhiStarGaussian', 'PhiStarGaussianComposition'};
measure_labels = {'P', 'C', '\Phi^{3}', '\phicomp', '\Phi*', '\Phi*comp', '\Phi*g', '\Phi*gcomp'};

%% Within-fly classification

results = cell(1, 2);

class_type = 'within';

for measure = 1 : length(measures)
    results{measure} = load([results_location measures{measure} '_svm_' class_type '.mat']);
end

% Average accuracies across flies
within = [];
for measure = 1 : length(results)
    nFlies = length(results{measure}.accuracy);
    
    accuracy = mean(results{measure}.accuracy);
    
    disp([measures{measure} ' - mean within-fly accuracy across N=13 flies: ' num2str(accuracy)]);
    
    within = cat(2, within, results{measure}.accuracy);
end

disp('====================================');

figure;
errorbar((1:length(results)), mean(within, 1), std(within, [], 1) / sqrt(size(within, 1)));
set(gca, 'XTick', (1:length(results)), 'XTickLabel', measure_labels);

figure;
plot((1:length(results)), within); legend(string(num2cell((1:13))), 'Location', 'southeast');
set(gca, 'XTick', (1:length(results)), 'XTickLabel', measure_labels);

%% Across-fly classification

results = cell(1, 2);

class_type = 'across';

for measure = 1 : length(measures)
    results{measure} = load([results_location measures{measure} '_svm_' class_type '.mat']);
end

across = zeros(size(results));
for measure = 1 : length(results)
    disp([measures{measure} ' - across-fly accuracy: ' num2str(results{measure}.accuracy)]);
    across(measure) = results{measure}.accuracy;
end

figure;
plot(across);
set(gca, 'XTick', (1:length(results)), 'XTickLabel', measure_labels);
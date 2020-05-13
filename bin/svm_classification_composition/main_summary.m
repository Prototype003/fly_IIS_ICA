%% Description

%{

Summarise SVM classification results

%}

%% Setup

results_location = 'results/';
measures = {'4ch_medianSplit_power', '4ch_medianSplit_coherence', '4ch_medianSplit_correlation', '2ch_phi3Composition_unpart', '3ch_phi3Composition_unpart', '4ch_phi3Composition_unpart', '4ch_phi3Value', '4ch_phiStarComposition', '4ch_phiStarGaussianComposition'};
measure_labels = {'P', 'C', 'r', '2ch\phi3comp', '3ch\phi3comp', '4ch\phi3comp', '4ch\Phi3', '\Phi*comp', '\Phi*gcomp'};

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

% Average accuracies across sets
%accuracy = permute(mean(results_mat, 1), [2 3 1]);

% Plot accuracy of each channel set, for each measure
figure; imagesc(accuracy, clim); c = colorbar; ylabel(c, 'accuracy');
set(gca, 'XTickLabel', measure_labels);
ylabel('channel set');

figure;
errorbar(mean(accuracy, 1), std(accuracy, [], 1), '.');
xlim([0 length(results)+1]); ylim([0.5 1]);
set(gca, 'XTick', (1: length(results)), 'XTickLabel', measure_labels);
ylabel('accuracy');

% Plot cumulative histogram for each measure
figure;
for measure = 1 : length(measures)
    cdfplot(mean(results{measure}.accuracy, 2)); % Average across subjects
    hold on;
end
legend(measure_labels, 'Location', 'southeast');

%% Stats

% Repeated measures ANOVA
t = array2table(accuracy);
rm = fitrm(t, ['accuracy1-accuracy' num2str(length(measures)) '~1']);
r = ranova(rm)
multcompare(rm, 'Time', 'ComparisonType', 'bonferroni')

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
figure; imagesc(accuracy, clim); c = colorbar; ylabel(c, 'accuracy');
set(gca, 'XTickLabel', measure_labels);
ylabel('channel set');

figure;
errorbar(mean(accuracy, 1), std(accuracy, [], 1), '.');
xlim([0 length(results)+1]); ylim([0.5 1]);
set(gca, 'XTick', (1: length(results)), 'XTickLabel', measure_labels);
ylabel('accuracy');

% Plot cumulative histogram for each measure
figure;
for measure = 1 : length(measures)
    cdfplot(mean(results{measure}.accuracy, 2)); % Average across subjects
    hold on;
end
legend(measure_labels, 'Location', 'southeast');

%% Stats

% Repeated measures ANOVA
t = array2table(accuracy);
rm = fitrm(t, ['accuracy1-accuracy' num2str(length(measures)) '~1']);
r = ranova(rm)
multcompare(rm, 'Time', 'ComparisonType', 'bonferroni')
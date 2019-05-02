%% Load accuracies

results_dir = 'results/';
results_file = ['phi3Composition_' constellation_type '_svm_' class_type.mat'];
load([results_dir results_file]);

%% Plot

% flies = (1:13);
% 
% order_labels = {'1st', '2nd', '3rd', '4th', '1/2/3/4th', '\Phi'};
% 
% % Plots average across flies, for all channel sets
% figure;
% imagesc(squeeze(mean(accuracies(flies, :, 1:6), 1)));
% c = colorbar; ylabel(c, 'accuracy');
% set(gca, 'XTickLabels', order_labels);
% xlabel('concept-orders used');
% ylabel('channel set');
% 
% % Plots average across channel sets, and average + std across flies
% figure;
% errorbar((1:6), squeeze(mean(mean(accuracies(flies, :, :), 2), 1)), squeeze(std(mean(accuracies(flies, :, :), 2), [], 1) / size(accuracies, 1)));
% %errorbar((1:6), squeeze(mean(mean(accuracies, 2), 1)), squeeze(std(mean(accuracies, 2), [], 1)));
% set(gca, 'XTickLabels', order_labels);
% xlabel('concept-orders used');
% ylabel('accuracy');
% xlim([0.5 6.5]);
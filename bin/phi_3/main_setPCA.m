%% Description

%{

Conducts PCA across channel sets (treats channel sets as dimensions)

Observations are big Phi for each channel set

Considerations:
    Observations - across flies, or PCA per fly?
    Observations - should include both wake/anest together
    Observations - include all trials? Or average before PCA?

%}

%% Load

load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM0.mat');
channel_sets = double(phis{1}.channel_sets);

%% Format for PCA

% Will include all flies, trials, conditions
flies = [1 2 3 4 5 6 7 8 9 10 12 13]; % Excluding 11


phi_values = [];

order_condition = zeros(numel(phis{1}.phis) / size(phis{1}.phis, 1) / size(phis{1}.phis, 2), 1);
order_fly = zeros(size(order_condition));
order_trial = zeros(size(order_condition));
row_counter = 1;
for condition = 1 : size(phis{1}.phis, 4)
    for fly = 1 : size(phis{1}.phis, 3)
        phi_values = cat(1, phi_values, permute(mean(phis{1}.phis(:, :, fly, condition), 2), [3 1 2 4]));
        order_condition(row_counter) = condition*2;
        order_fly(row_counter) = fly;
        order_trial(row_counter) = trial;
        row_counter = row_counter + 1;
%         for trial = 1 : size(phis{1}.phis, 2)
%             phi_values = cat(1, phi_values, permute(phis{1}.phis(:, trial, fly, condition), [2 1 3 4]));
%             order_condition(row_counter) = condition*2;
%             order_fly(row_counter) = fly;
%             order_trial(row_counter) = trial;
%             row_counter = row_counter + 1;
%         end
    end
end

% PCA
[coeff, score, latent, tsquared, explained] = pca(zscore(phi_values));
size(score)

%mapped = tsne(zscore(phi_values));

% Plot

figure; colormap('copper');
% scatter(score(1:size(score, 1)/2, 1), score(1:size(score, 1)/2, 2), [], (1:13), 'o'); hold on;
% scatter(score(1+size(score, 1)/2:end, 1), score(1+size(score, 1)/2:end, 2), [], (1:13), 'x');
scatter(score(:, 1), score(:, 2), 100, order_condition(find(order_condition)), '.');
%scatter3(score(:, 1), score(:, 2), score(:, 3), [], order_condition(find(order_condition)), '.');
%scatter(mapped(:, 1), mapped(:, 2), [], order_condition, '.');
% load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar.mat')
% phis1 = phis;
% load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar_phiToolbox.mat')
% phis2 = phis;


%% Plot correlations

setSize = 3;
plot_metrics = {...
    'hs',...
    'hconds',...
    'mis',...
    'mi_stars',...
    'phi_stars',...
    'phi_stars_normalised'};
plot_titles = {...
    'H',...
    'Hcond',...
    'MI',...
    'MI*',...
    'Phi*'...
    'Phi* norm'};
plot_positions = [1 4 2 5 3 6];
subplot_dimensions = [2 3];

figure;
for plot_counter = 1 : length(plot_metrics)
    subplot(subplot_dimensions(1), subplot_dimensions(2), plot_positions(plot_counter));
    scatter(phis1{setSize}.(plot_metrics{plot_counter})(:), phis2{setSize}.(plot_metrics{plot_counter})(:), 0.001, '.', 'MarkerEdgeAlpha', 0.01);
    axis tight
    title(plot_titles{plot_counter});
    if plot_counter == 1
        xlabel('Feb2014'); ylabel('phitoolbox');
    end
end

%% Find set where phi-star difference is maximal

setSize = 2;

phi_diffs = abs(phis1{setSize}.phi_stars - phis2{setSize}.phi_stars);
mi_diffs = abs(phis1{setSize}.mis - phis2{setSize}.mis);

[diff_max, diff_position] = max(phi_diffs(:));
[set_counter, trial, fly, condition, tau] = ind2sub(size(phi_diffs), diff_position);

% [diff_max, diff_position] = max(mi_diffs(:));
% [set_counter, trial, fly, condition, tau] = ind2sub(size(phi_diffs), diff_position);

set = phis{setSize}.channel_sets(set_counter, :);

figure;
% Show data point with maximal difference (on phi-star plot)
subplot(1, 2, 1);
scatter(phis1{setSize}.mis(:), phis2{setSize}.mis(:), '.');
hold on;
scatter(phis1{setSize}.mis(set_counter, trial, fly, condition, tau), phis2{setSize}.mis(set_counter, trial, fly, condition, tau), 'r.')
axis tight
title('MI');

% Show data point with maximal difference (on phi-star plot)
subplot(1, 2, 2);
scatter(phis1{setSize}.phi_stars(:), phis2{setSize}.phi_stars(:), '.');
hold on;
scatter(phis1{setSize}.phi_stars(set_counter, trial, fly, condition, tau), phis2{setSize}.phi_stars(set_counter, trial, fly, condition, tau), 'r.')
axis tight
title('Phi*');


%% Get fly data corresponding to max difference

fly_data = load('workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat');
fly_data = fly_data.fly_data;

set_data = fly_data(:, set, trial, fly, condition)';

%% Calculate phi values

taus = [4 8 16];

% Covariances
[cov_present_present, cov_present_past, cov_past_past] = Cov_comp_sample(set_data, taus(tau));
[cpa, cpa_pr, cpr] = Cov_comp(set_data, taus(tau));

% phi_toolbox_Feb2014
%phi_Feb2014 = phistar_mip2(cpa, cpa_pr', cpr, set);

% phi_toolbox
phi_phiToolbox = phistar_mip(cpa, cpa_pr, cpr, set);

%%

channel_data = struct();
channel_data.phi_diff = set_data;
channel_data.phi_diff_tau = taus(tau);

channel_data.mi_diff = set_data;
channel_data.mi_diff_tau = taus(tau);

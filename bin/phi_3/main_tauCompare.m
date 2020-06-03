%% Description

%{
Plots big-phi at each tau level
%}

%% Setup

%load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials_nComponents4_phithree_nChannels4_globalTPM0.mat');
load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials_nComponents4_phithree_nChannels4_globalTPM0_binOffset1.mat');

%% Plot big-phi
% 1 plot per fly, comparing all trials of wake vs anesthesia

cond_colours = {'r', 'b'};

big_phis = phis{1}.phis;

figure;
for fly = 1 : size(big_phis, 3)
    subplot(4, 4, fly);
    for condition = 1 : size(big_phis, 4)
        % Plot individual trials
        plot(phis{1}.taus, squeeze(big_phis(:, :, fly, condition, :)), [cond_colours{condition} '-.']);
        hold on;
        % Plot average across trials
        plot(phis{1}.taus, squeeze(mean(big_phis(:, :, fly, condition, :), 2)), cond_colours{condition}, 'LineWidth', 3);
    end
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('\tau (ms)');
    ylabel('\Phi');
    title(['fly' num2str(fly)]);
    
    %ylim([min(big_phis(:)) max(big_phis(:))]);
end

% Population (all flies)
subplot(4, 4, fly+1);
for condition = 1 : size(big_phis, 4)
    % Plot individual flies
    plot(phis{1}.taus, squeeze(mean(big_phis(:, :, :, condition, :), 2)), [cond_colours{condition} '-.']);
    hold on;
    % Plot average across flies
    plot(phis{1}.taus, squeeze(mean(mean(big_phis(:, :, :, condition, :), 2), 3)), cond_colours{condition}, 'LineWidth', 3);
end
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('\tau (ms)');
ylabel('\Phi');
title(['all flies']);

%% Plot big-phi (wake - anest)

figure;
for fly = 1 : size(big_phis, 3)
    subplot(4, 4, fly);
    data = big_phis(:, :, fly, 1, :) ./ big_phis(:, :, fly, 2, :);
    % Plot individual trials
    plot(phis{1}.taus, squeeze(data), 'k-.');
    hold on;
    % Plot average across trials
    plot(phis{1}.taus, squeeze(mean(data, 2)), 'k', 'LineWidth', 3);
    
    set(gca, 'XScale', 'log');
    xlabel('tau (ms)');
    ylabel('\Delta\Phi');
    title(['fly' num2str(fly)]);
end

% Population (all flies)
subplot(4, 4, fly+1)
data = big_phis(:, :, :, 1, :) ./ big_phis(:, :, :, 2, :);
% Plot individuals flies
plot(phis{1}.taus, squeeze(mean(data, 2)), 'k-.');
hold on;
% Plot average across flies
plot(phis{1}.taus, squeeze(mean(mean(data, 2), 3)), 'k', 'LineWidth', 3);
set(gca, 'XScale', 'log');
xlabel('\tau (ms)');
ylabel('\Delta\Phi');
title(['all flies']);
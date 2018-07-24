%% DESCRIPTION

%{
Paper figure 2

Script for figure 1 is tpm_worked_example.m

Original script: main_phiplots.m
%}

%% SETUP

global_tpm = 1;

%% SETUP

flies = (1:13); %[1 2 3 4 5 6 7 9 10 11 12 13];

fontsize = 11;

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

filter_percent = 5; % 0 to 100 %

tau_colours = {'r', 'g', 'b'};
tau_alphas = [1 0.8 0.5];

bin_location = '../';
addpath(bin_location);

%% LOAD

[phi_threes, measure_strings{1}] = phi_load('phi_three', global_tpm, bin_location);

[phi_stars, measure_strings{2}] = phi_load('phi_star', 0, bin_location);

% %% Rename some stuff
% 
% for nChannels_counter = 1 : length(phi_threes)
%     phi_threes{nChannels_counter}.phis = phi_threes{nChannels_counter}.phi_threes;
%     phi_stars{nChannels_counter}.phis = phi_stars{nChannels_counter}.phi_stars;
% end

%% Reset phi variable to original state
% for nChannels_counter = 1 : numel(phi_threes)
%     phi_threes{nChannels_counter}.phis = phi_threes{nChannels_counter}.phis_raw;
%     phi_stars{nChannels_counter}.phis = phi_stars{nChannels_counter}.phis_raw;
% end

%% Average across trials and flies, calculate deltas and associated standard error

% Stats
for nChannels_counter = 1 : numel(phi_threes)
    % Raw
    phi_threes{nChannels_counter}.phis_raw = phi_threes{nChannels_counter}.phis(:, :, flies, :, :);
    phi_stars{nChannels_counter}.phis_raw = phi_stars{nChannels_counter}.phis(:, :, flies, :, :);
    
    % Average 
    phi_threes{nChannels_counter}.phis = squeeze(mean(mean(phi_threes{nChannels_counter}.phis(:, :, flies, :, :), 2), 3));
    phi_stars{nChannels_counter}.phis = squeeze(mean(mean(phi_stars{nChannels_counter}.phis(:, :, flies, :, :), 2), 3));
    
    % Standard error (across flies)
    phi_threes{nChannels_counter}.phis_std = squeeze(std(mean(phi_threes{nChannels_counter}.phis_raw(:, :, flies, :, :), 2), [], 3)) / sqrt(length(flies));
    phi_stars{nChannels_counter}.phis_std = squeeze(std(mean(phi_stars{nChannels_counter}.phis_raw(:, :, flies, :, :), 2), [], 3)) / sqrt(length(flies));
    
    % Delta (air - iso) (absolute)
    phi_threes{nChannels_counter}.phis_delta = squeeze(mean(mean(phi_threes{nChannels_counter}.phis_raw(:, :, flies, 1, :) - phi_threes{nChannels_counter}.phis_raw(:, :, flies, 2, :), 2), 3));
    phi_stars{nChannels_counter}.phis_delta = squeeze(mean(mean(phi_stars{nChannels_counter}.phis_raw(:, :, flies, 1, :) - phi_stars{nChannels_counter}.phis_raw(:, :, flies, 2, :), 2), 3));
    
    % Delta std (air - iso)
    phi_threes{nChannels_counter}.phis_delta_std = squeeze(std(mean(phi_threes{nChannels_counter}.phis_raw(:, :, flies, 1, :) - phi_threes{nChannels_counter}.phis_raw(:, :, flies, 2, :), 2), [], 3)) / sqrt(length(flies));
    phi_stars{nChannels_counter}.phis_delta_std = squeeze(std(mean(phi_stars{nChannels_counter}.phis_raw(:, :, flies, 1, :) - phi_stars{nChannels_counter}.phis_raw(:, :, flies, 2, :), 2), [], 3)) / sqrt(length(flies));
    
%     % Delta (air - iso / air) (relative)
%     phi_threes{nChannels_counter}.phis_delta = squeeze(mean(mean((phi_threes{nChannels_counter}.phi_threes(:, :, flies, 1, :) - phi_threes{nChannels_counter}.phi_threes(:, :, flies, 2, :)) ./ phi_threes{nChannels_counter}.phi_threes(:, :, flies, 1, :), 2), 3));
%     phi_stars{nChannels_counter}.phis_delta = squeeze(mean(mean((phi_stars{nChannels_counter}.(star_metric)(:, :, flies, 1, :) - phi_stars{nChannels_counter}.(star_metric)(:, :, flies, 2, :)) ./ phi_stars{nChannels_counter}.(star_metric)(:, :, flies, 1, :), 2), 3));
%     
%     % Delta std (air - iso)
%     phi_threes{nChannels_counter}.phis_delta_std = squeeze(std(mean((phi_threes{nChannels_counter}.phi_threes(:, :, flies, 1, :) - phi_threes{nChannels_counter}.phi_threes(:, :, flies, 2, :)) ./phi_threes{nChannels_counter}.phi_threes(:, :, flies, 1, :), 2), [], 3)) / sqrt(length(flies));
%     phi_stars{nChannels_counter}.phis_delta_std = squeeze(std(mean((phi_stars{nChannels_counter}.(star_metric)(:, :, flies, 1, :) - phi_stars{nChannels_counter}.(star_metric)(:, :, flies, 2, :)) ./ phi_stars{nChannels_counter}.(star_metric)(:, :, flies, 1, :), 2), [], 3)) / sqrt(length(flies));
end

%% Average across trials and channel sets, at each set size
% channel_ticks holds the set value at which nChannels increments

% nChannels x taus x conditions x flies
phis_averaged = struct();
phis_averaged.threes = zeros(length(phi_threes), size(phi_threes{1}.phis_raw, 5), size(phi_threes{1}.phis_raw, 4), size(phi_threes{1}.phis_raw, 3));
phis_averaged.stars = zeros(length(phi_threes), size(phi_threes{1}.phis_raw, 5), size(phi_threes{1}.phis_raw, 4), size(phi_threes{1}.phis_raw, 3));

for tau = 1 : size(phi_threes{1}.phis_raw, 5)
    for nChannels = 1 : length(phi_threes)
        for condition = 1 : size(phi_threes{1}.phis_raw, 4)
            phis_averaged.threes(nChannels, :, :, :) = permute(mean(mean(phi_threes{nChannels}.phis_raw, 2), 1), [5 4 3 2 1]);
            phis_averaged.stars(nChannels, :, :, :) = permute(mean(mean(phi_stars{nChannels}.phis_raw, 2), 1), [5 4 3 2 1]);
        end
    end
end

%% Significance tests on averages

stat_metric = 'threes'; % Or 'stars'

%% Post-hoc tests on network size, using averages
foo = squeeze(mean(mean(phis_averaged.(stat_metric), 2), 3));
[p, h, stats] = signrank(foo(2,:), foo(3,:), 'tail', 'left', 'method', 'exact');

%% Post-hoc tests on lag, using averages
tmp = squeeze(mean(phis_averaged.(stat_metric), 3));
foo = zeros(size(phis_averaged.(stat_metric), 2), size(phis_averaged.(stat_metric), 4));
for tau_counter = 1 : size(phis_averaged.(stat_metric), 1)
    for fly = 1 : size(phis_averaged.(stat_metric), 4)
        % Weighted average
        foo(tau_counter, fly) =...
            ((tmp(1, tau_counter, fly) * 105)...
            +(tmp(2, tau_counter, fly) * 455)...
            +(tmp(3, tau_counter, fly) * 1365))...
            /(105+455+1365);
    end
end
[p, h, stats] = signrank(foo(2,:), foo(3,:), 'method', 'exact');

%% Plot

figure;
set(gcf, 'Position', [0 0 2100/3 700]);
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
% set(gcf, 'RendererMode', 'manual');
% set(gcf, 'Renderer', 'painters');

% Plot marker shapes/lines/colours
setSize_caps = [2 4 6];
setSize_widths = [0.25 0.5 1];
setSize_markers = 'x^o';
setSize_lines = {':', '--', '-'};
setSize_offsets = [-0.15 0 0.15];
taus = [4 8 16];

plot_metrics = {'threes', 'stars'};
metric_labels = {'^{3.0}', '*'};
metric_lims = struct();
metric_limits.threes = [0.03 0.03 0.015];
metric_limits.stars = [0.005 0.005 0.002]; % 0.005 0.005 0.002 % For Gaussian assumption values
metric_limits.stars = [0.05 0.05 0.002];
metric_exponents.threes = -2;
metric_exponents.stars = -3; % For Gaussian assumption values
metric_exponents.stars = -2;

yticks = struct();
yticks.threes = [0 metric_limits.threes/2 metric_limits.threes];
yticks.stars = [0 metric_limits.stars/2 metric_limits.stars];


% subplot positioning
ySpacing = 0.1;
xSpacing = 0.2;

yPortion = 0.7;
xPortion = 0.7;

textbox_width = 0.03;
text_labels = 'adbecf';

heights = [1/3 1/3 1/3] * yPortion; % y-axis doubles each time, x + 2x + 4x = 1
widths = [3/5 3/5] * (xPortion - xSpacing);

yStarts = (1-yPortion)/2 + fliplr([0 cumsum(fliplr(heights(2:end)))]);
xStarts = (1-xPortion)/3 + (0 : widths(1) : 1);
xStarts(2) = xStarts(2) + xSpacing;

for metric_counter = 1 : length(plot_metrics)
    plot_metric = plot_metrics{metric_counter};
    
    % Air phi values
    subplot(3, length(plot_metrics), metric_counter); % Each column repeats the same plots, for a different measure
    for setSize_counter = 1 : size(phis_averaged.(plot_metric), 1)
        values = mean(phis_averaged.(plot_metric)(setSize_counter, :, 1, :), 4);
        values_std = std(phis_averaged.(plot_metric)(setSize_counter, :, 1, :), [], 4) / sqrt(size(phis_averaged.(plot_metric), 4));
        errorbar((1:length(taus)) + setSize_offsets(setSize_counter),...
            values,...
            values_std,...
            ['k' setSize_markers(setSize_counter) setSize_lines{setSize_counter}],...
            'CapSize', setSize_caps(setSize_counter),...
            'LineWidth', setSize_widths(setSize_counter),...
            'MarkerFaceColor', 'k',...
            'MarkerSize', 4);
        hold on;
    end
    %legend('    2ch', '    3ch', '    4ch', 'Location', 'southeast');
    set(gca, 'Position', [xStarts(metric_counter) yStarts(1)+ySpacing, widths(metric_counter), heights(1)-ySpacing]);
    axis_defaults(gca);
    xlim([0.75 3.25]);
    ylim([0 metric_limits.(plot_metric)(1)]);
    max_y = ylim; max_y = max_y(end);
    set(gca, 'YTick', [0 max_y/2 max_y]);
    ax = gca;
    ax.YAxis.Exponent = metric_exponents.(plot_metric);
    set(gca, 'XTick', (1:length(taus)));
    set(gca, 'XTickLabel', taus);
    y = ylabel(['\Phi' metric_labels{metric_counter}], 'rotation', 90);
    title(['  \Phi' metric_labels{metric_counter} '\newlinewake']);
    ax_pos = get(gca, 'Position');
    axes('Visible', 'off', 'Position', [ax_pos(1)-0.05 ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
    text(0, 0.2, text_labels(metric_counter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');
    
    
    % Iso phi values
    subplot(3, length(plot_metrics), metric_counter + 2);
    for setSize_counter = 1 : size(phis_averaged.(plot_metric), 1)
        values = mean(phis_averaged.(plot_metric)(setSize_counter, :, 2, :), 4);
        values_std = std(phis_averaged.(plot_metric)(setSize_counter, :, 2, :), [], 4) / sqrt(size(phis_averaged.(plot_metric), 4));
        errorbar((1:length(taus)) + setSize_offsets(setSize_counter),...
            values,...
            values_std,...
            ['k' setSize_markers(setSize_counter) setSize_lines{setSize_counter}],...
            'CapSize', setSize_caps(setSize_counter),...
            'LineWidth', setSize_widths(setSize_counter),...
            'MarkerFaceColor', 'k',...
            'MarkerSize', 4);
        hold on;
    end
    set(gca, 'Position', [xStarts(metric_counter) yStarts(2)+ySpacing, widths(metric_counter), heights(2)-ySpacing]);
    axis_defaults(gca);
    xlim([0.75 3.25]);
    ylim([0 metric_limits.(plot_metric)(1)]);
    max_y = ylim; max_y = max_y(end);
    set(gca, 'YTick', [0 max_y/2 max_y]);
    ax = gca;
    ax.YAxis.Exponent = metric_exponents.(plot_metric);
    set(gca, 'XTick', (1:length(taus)));
    set(gca, 'XTickLabel', taus);
    y = ylabel(['\Phi' metric_labels{metric_counter}], 'rotation', 90);
    title('anes');
    ax_pos = get(gca, 'Position');
    axes('Visible', 'off', 'Position', [ax_pos(1)-0.05 ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
    text(0, 0.2, text_labels(metric_counter+2), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');
    
    % Air - Iso phi values
    subplot(3, length(plot_metrics), metric_counter + 4);
    for setSize_counter = 1 : size(phis_averaged.(plot_metric), 1)
        values = mean(phis_averaged.(plot_metric)(setSize_counter, :, 1, :) - phis_averaged.(plot_metric)(setSize_counter, :, 2, :), 4);
        values_std = std(phis_averaged.(plot_metric)(setSize_counter, :, 1, :) - phis_averaged.(plot_metric)(setSize_counter, :, 2, :), [], 4) / sqrt(size(phis_averaged.(plot_metric), 4));
        errorbar((1:length(taus)) + setSize_offsets(setSize_counter),...
            values,...
            values_std,...
            ['k' setSize_markers(setSize_counter) setSize_lines{setSize_counter}],...
            'CapSize', setSize_caps(setSize_counter),...
            'LineWidth', setSize_widths(setSize_counter),...
            'MarkerFaceColor', 'k',...
            'MarkerSize', 4);
        hold on;
    end
    set(gca, 'Position', [xStarts(metric_counter) yStarts(3)+ySpacing, widths(metric_counter), heights(3)-ySpacing]);
    axis_defaults(gca);
    mins.threes = -0.003;
    mins.stars = 0;
    xlim([0.75 3.25]);
    ylim([mins.(plot_metric) metric_limits.(plot_metric)(3)]);
    max_y = ylim; max_y = max_y(end);
    set(gca, 'YTick', [0 max_y/2 max_y]);
    ax = gca;
    ax.YAxis.Exponent = metric_exponents.(plot_metric);
    set(gca, 'XTick', (1:length(taus)));
    set(gca, 'XTickLabel', taus);
    y = ylabel(['\Delta\Phi' metric_labels{metric_counter}], 'rotation', 90);
    title('wake - anes');
    xlabel('\tau (ms)');
    ax_pos = get(gca, 'Position');
    axes('Visible', 'off', 'Position', [ax_pos(1)-0.05 ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
    text(0, 0.2, text_labels(metric_counter+4), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');
    
end

%% Print figure

% figure_name = 'fig2';
% 
% print(figure_name, '-dsvg'); % SVG
% print(figure_name, '-dpdf', '-bestfit'); % PDF
% print(figure_name, '-dpng'); % PNG
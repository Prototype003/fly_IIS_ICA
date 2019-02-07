%% DESCRIPTION

%{

Figure 5 - correlation between mean 2ch mean values with 3/4ch values
2 x 4 (magnitude/class by [3/4ch by within/across])

%}

%% Setup

measure = 'phi_three'; % 'phi_three' or 'phi_star' or 'phi_star_gaussian' or 'phi_SI'
tau = 1; % 1 = 4ms; 2 = 8ms; 3 = 16ms
if tau == 1
    tau_string = '4';
elseif tau == 2
    tau_string = '8';
elseif tau == 3
    tau_string = '16';
end

freq_range = (1:42); %(1:83); % corresponding to ~5Hz and ~10Hz, check the 'frequencies' vector
freq_range_string = '0-5Hz'; %'0-10Hz';

fontsize = 11; % Used for drawing label letters

bin_location = '../';
addpath(bin_location);

results_directory = [bin_location 'workspace_results/'];

%% Load phi values (WITHIN)

% We'll not rename this because the filesize is big and previous figures use this as well
[phis, measure_string] = phi_load(measure, 0, bin_location);

%% Load accuracy results for phi

results_directory = '../workspace_results/';

if strcmp(measure, 'phi_three')
    % Load ACROSS
    results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification_across1.mat';
    load([results_directory results_filename]);
    % Fix python indexing
    for nChannels_counter = 1 : length(accuracies)
        accuracies{nChannels_counter}.channel_sets = accuracies{nChannels_counter}.channel_sets + 1;
    end
    accuracies_a = accuracies;
    % Load WITHIN
    results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
    load([results_directory results_filename]);
    % Fix python indexing
    python_indexing = [1 1 0];
    for nChannels_counter = 1 : length(accuracies)
        accuracies{nChannels_counter}.channel_sets = accuracies{nChannels_counter}.channel_sets + python_indexing(nChannels_counter);
    end
    accuracies_w = accuracies;
else % strcmp(measure, 'phi_star')
    % Load ACROSS
    results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification_across1.mat';
    load([results_directory results_filename]);
    accuracies_a = accuracies;
    % Load WITHIN
    results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
    load([results_directory results_filename]);
    accuracies_w = accuracies;
end

%% Get mean values per channel

channels = (1 : max(accuracies_w{1}.channel_sets(:)));
channel_means = cell(length(phis), 1); % We will save channel means so we can use them for selecting large sets of channels

% Storage of channel means
for nChannels_counter = 1 : length(accuracies)
    channel_sets = accuracies_w{nChannels_counter}.channel_sets;
    
    % Phi values
    tmp_size = size(phis{nChannels_counter}.phis(:, :, :, :, tau));
    tmp_size(1) = length(channels);
    phis{nChannels_counter}.channel_sums = zeros(tmp_size);
    phis{nChannels_counter}.channel_means = zeros(tmp_size);
    phis{nChannels_counter}.set_counters = zeros(tmp_size);
    
    % Accuracies
    tmp_size = size(accuracies_w{nChannels_counter}.accuracies);
    tmp_size(1) = length(channels);
    accuracies_w{nChannels_counter}.channel_sums = zeros(tmp_size);
    accuracies_w{nChannels_counter}.channel_means = zeros(tmp_size);
    accuracies_w{nChannels_counter}.set_counters = zeros(tmp_size);
    accuracies_a{nChannels_counter}.channel_sums = zeros(tmp_size(1), 1);
    accuracies_a{nChannels_counter}.channel_means = zeros(tmp_size(1), 1);
    accuracies_a{nChannels_counter}.set_counters = zeros(tmp_size(1), 1);
    
    % Sum values for each channel (sum across networks which contain the channel)
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                
                phis{nChannels_counter}.channel_sums(channel, :, :, :) = phis{nChannels_counter}.channel_sums(channel, :, :, :) + phis{nChannels_counter}.phis(channel_set, :, :, :, tau);
                phis{nChannels_counter}.set_counters(channel, :, :, :) = phis{nChannels_counter}.set_counters(channel, :, :, :) + 1;
                
                accuracies_w{nChannels_counter}.channel_sums(channel, :) = accuracies_w{nChannels_counter}.channel_sums(channel, :) + accuracies_w{nChannels_counter}.accuracies(channel_set, :);
                accuracies_w{nChannels_counter}.set_counters(channel, :) = accuracies_w{nChannels_counter}.set_counters(channel, :) + 1;
                
                accuracies_a{nChannels_counter}.channel_sums(channel, :) = accuracies_a{nChannels_counter}.channel_sums(channel, :) + accuracies_a{nChannels_counter}.accuracies(channel_set, :);
                accuracies_a{nChannels_counter}.set_counters(channel, :) = accuracies_a{nChannels_counter}.set_counters(channel, :) + 1;
                
            end
        end
    end
    
    % Average values
    phis{nChannels_counter}.channel_means = phis{nChannels_counter}.channel_sums ./ phis{nChannels_counter}.set_counters;
    accuracies_w{nChannels_counter}.channel_means = accuracies_w{nChannels_counter}.channel_sums ./ accuracies_w{nChannels_counter}.set_counters;
    accuracies_a{nChannels_counter}.channel_means = accuracies_a{nChannels_counter}.channel_sums ./ accuracies_a{nChannels_counter}.set_counters;
    
end

%% Make figure

figure;
set(gcf, 'Position', [0 0 2100/4 600]);
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing

titles = {'3ch', '4ch'};
class_titles = {'within', 'across'};

textbox_width = 0.03;
text_labels = 'abcdefgh';

ySpacing = 0.1;
xSpacing = 0.05;

yPortion = 0.7;
xPortion = 0.7;

heights = [1 1 1] * yPortion ./ 3;
widths = ([1 1] * xPortion ./ 2);

y_offset = 0;
yStarts = linspace(0, 1, length(heights)+1);
yStarts = fliplr(yStarts(2:end) - heights(1));
yStarts = yStarts + y_offset;

x_offset = 0.01;
xStarts = linspace(0, 1, length(widths)+1);
xStarts = xStarts(2:end) - widths(1);
xStarts = xStarts + x_offset;

xCounter = 1;
yCounter = 1;
subplot_counter = 1;

colours = 'rb';
shapes = 'ox'; shapes = '.x';

xlims = [...
    0.001 0.004
    0.001 0.004
    56.5 60.5;
    56.5 60.5;
    54 63;
    54 63];

% Phi-value correlation
for nChannels_counter = 2 : length(phis)
    subplot(size(xlims, 1) / 2, length(phis) - 1, subplot_counter);
    set(gca, 'Position', [xStarts(xCounter) yStarts(yCounter), widths(xCounter)-xSpacing, heights(yCounter)-ySpacing]);
    
    channel_sets = phis{nChannels_counter}.channel_sets;
    
    values = phis{nChannels_counter}.phis(:, :, :, :, tau);
    channel_means = phis{1}.channel_means(:, :, :, :);
    
    % Build values from average of single channel average values
    values_rebuilt = zeros(size(values));
    for channel_set = 1 : size(values, 1)
        channels = channel_sets(channel_set, :);
        values_rebuilt(channel_set, :, :, :) = mean(channel_means(channels, :, :, :), 1);
    end
    
    % Store
    phis{nChannels_counter}.phis_rebuilt = values_rebuilt;
    
    % Plots wake and anest
    for condition = 1 : size(values_rebuilt, 4)
        flies = (1:13);
        %flies = flies(flies ~= 11);
        x = mean(mean(values_rebuilt(:, :, flies, condition), 2), 3); % Average across trials and flies
        y = mean(mean(values(:, :, flies, condition), 2), 3);
        scatter(x(:), y(:), [], colours(condition), shapes(condition)); hold on;
    end
    
    % Plots wake - anest
    %     flies = [1 2 3 4 5 6 7 8 9 10 12 13];
    %     x = mean(values_rebuilt(:, :, flies, 1) - values_rebuilt(:, :, flies, 2), 2);
    %     y = mean(values(:, :, flies, 1) - values(:, :, flies, 2), 2);
    %     scatter(x(:), y(:), [], colours(condition), shapes(condition)); hold on;
    
    title([num2str(nChannels_counter+1) 'ch']);
    
    if nChannels_counter == 2
        ylabel(measure_string);
        xlabel(['mean 2ch ' measure_string]);
    end
    
    handle = gca;
    handle.XRuler.Exponent = 0;
    
    axis_defaults(gca);
    
    set(gca, 'Box', 'on');
    tmp = xlim;
    tmp_range = tmp(2) - tmp(1);
    tmp_adjust = tmp_range * 0.15;
    %set(gca, 'XTick', linspace(tmp(1)+tmp_adjust, tmp(2)-tmp_adjust, 3));
    
    %axis equal;
    xlim(xlims(subplot_counter, :));
    
    % Add letter label to plot
    ax_pos = get(gca, 'Position');
    axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
    text(0, 0.2, text_labels(subplot_counter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');
    
    subplot_counter = subplot_counter + 1;
    xCounter = xCounter + 1;
    
end
yCounter = yCounter + 1;

% Classification correlation
for class_type_counter = 1 : length(class_titles)
    class_type = class_titles{class_type_counter};
    xCounter = 1;
    
    for nChannels_counter = 2 : length(phis)
        subplot(size(xlims, 1) / 2, length(phis) - 1, subplot_counter);
        set(gca, 'Position', [xStarts(xCounter) yStarts(yCounter), widths(xCounter)-xSpacing, heights(yCounter)-ySpacing]);
        
        channel_sets = phis{nChannels_counter}.channel_sets;
        
        if strcmp(class_type, 'within')
            values = accuracies_w{nChannels_counter}.accuracies;
            channel_means = accuracies_w{1}.channel_means;
        else % strcmp(class_type, 'across');
            values = accuracies_a{nChannels_counter}.accuracies;
            channel_means = accuracies_a{1}.channel_means;
        end
        
        % Build values from average of single channel average values
        values_rebuilt = zeros(size(values));
        for channel_set = 1 : size(channel_sets, 1)
            channels = channel_sets(channel_set, :);
            values_rebuilt(channel_set, :) = mean(channel_means(channels, :), 1);
        end
        
        % Store
        if strcmp(class_type, 'within')
            accuracies_w{nChannels_counter}.values_rebuilt = values_rebuilt;
        else % strcmp(class_type, 'across');
            accuracies_a{nChannels_counter}.values_rebuilt = values_rebuilt;
        end
        
        % Plot
        x = mean(values_rebuilt, 2);
        y = mean(values, 2);
        scatter(x(:), y(:), [], 'k.'); hold on;
        
        if nChannels_counter == 2
            ylabel('class. acc. %');
            xlabel('mean 2ch acc. %');
        end
        
        axis_defaults(gca);
        
        set(gca, 'Box', 'on');
        
        %axis equal;
        xlim(xlims(subplot_counter, :));
        axis tight;
        
        % Add letter label to plot
        ax_pos = get(gca, 'Position');
        axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
        text(0, 0.2, text_labels(subplot_counter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');
        
        % Add row heading
        if nChannels_counter == 3
            axes('Visible', 'off', 'Position', [(xStarts(xCounter-1)+widths(xCounter-1)+xSpacing) +yStarts(yCounter)+heights(yCounter)-0.05, 0.1 0.1]); axis_defaults(gca);
            text(0, 0, class_type, 'FontSize', fontsize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
        
        subplot_counter = subplot_counter + 1;
        xCounter = xCounter + 1;
    end
    yCounter = yCounter + 1;
end

%% Regress 3,4ch big-phi onto 2ch big-phi

% Will ignore effect of condition - group conditions together

% Make tables
tables = cell(size(phis));
for nChannels_counter = 2 : length(phis)
    values = double(permute(mean(phis{nChannels_counter}.phis(:, :, :, :, tau), 2), [1 3 4 2]));
    values_rebuilt = double(permute(mean(phis{nChannels_counter}.phis_rebuilt, 2), [1 3 4 2]));
    set_ids = zeros(size(values));
    fly_ids = zeros(size(values));
    conditions = zeros(size(values));
    
    for set_id = 1 : size(values, 1)
        set_ids(set_id, :, :) = set_id;
    end
    for fly = 1 : size(values, 2)
        fly_ids(:, fly, :) = fly;
    end
    for condition = 1 : size(values, 3)
        conditions(:, :, condition) = condition;
    end
    
    tables{nChannels_counter} = table(...
        values(:), values_rebuilt(:), set_ids(:), fly_ids(:), conditions(:),...
        'VariableNames', {'value', 'value_rebuilt', 'set_id', 'fly_id', 'condition'});
    
    tables{nChannels_counter}.fly_id = nominal(tables{nChannels_counter}.fly_id);
    tables{nChannels_counter}.set_id = nominal(tables{nChannels_counter}.set_id);
    tables{nChannels_counter}.condition = nominal(tables{nChannels_counter}.condition);
end

nChannels_counter = 3;

% Build model
model_spec = 'value ~ value_rebuilt + (1|fly_id) + (1|fly_id:set_id)';
model_full = fitlme(tables{nChannels_counter}, model_spec);

model_null_spec = 'value ~ 1 + (1|fly_id) + (1|fly_id:set_id)';
model_null = fitlme(tables{nChannels_counter}, model_null_spec);

compare(model_null, model_full)

%% Plot

%figure;

nChannels_counter = 3;

channel_sets = phis{nChannels_counter}.channel_sets;

values = phis{nChannels_counter}.phis(:, :, :, :, tau);
channel_means = phis{1}.channel_means(:, :, :, :);

% Build values from average of single channel average values
values_rebuilt = zeros(size(values));
for channel_set = 1 : size(values, 1)
    channels = channel_sets(channel_set, :);
    values_rebuilt(channel_set, :, :, :) = mean(channel_means(channels, :, :, :), 1);
end

% Store
phis{nChannels_counter}.phis_rebuilt = values_rebuilt;

% Plots wake and anest
xs = []; ys = [];
legend_items = [];
for condition = 1 : size(values_rebuilt, 4)
    flies = (1:13);
    %flies = flies(flies ~= 11);
    x = mean(mean(values_rebuilt(:, :, flies, condition), 2), 3); % Average across trials and flies
    y = mean(mean(values(:, :, flies, condition), 2), 3);
    legend_items(condition) = scatter(x(:), y(:), [], colours(condition), '.'); hold on;
    xs = [xs; x(:)]; ys = [ys; y(:)];
    %Plot line of best fit
    p = polyfit(x, y, 1);
    best_fit = polyval(p, x);
    plot(x, best_fit, 'k', 'LineWidth', 2);
    p
end
legend(legend_items, 'wake', 'anest', 'Location', 'northwest');

handle = gca;
handle.XRuler.Exponent = 0;
set(gca, 'XTick', [0.001 0.0025 0.004], 'YTick', [0.015 0.035 0.055]);
%axis([0.0009 0.005 -0.002 0.06]);
axis([floor(1000*min(xs))/1000 ceil(1000*max(xs))/1000 floor(1000*min(ys))/1000 ceil(1000*max(ys))/1000])
axis([0.001 0.004 0.013 0.06]);
axis square

% Plots wake - anest
% flies = [1 2 3 4 5 6 7 8 9 10 11 12 13];
% x = mean(mean(values_rebuilt(:, :, flies, 1) - values_rebuilt(:, :, flies, 2), 2), 3);
% y = mean(mean(values(:, :, flies, 1) - values(:, :, flies, 2), 2), 3);
% scatter(x(:), y(:), [], colours(condition), shapes(condition)); hold on;
% xs = x(:); ys = y(:);
% 
% % Plot line of best fit
% p = polyfit(xs, ys, 1);
% best_fit = polyval(p, xs);
% plot(xs, best_fit);

%% Regress 3,4ch classification onto 2ch classification (within)

% Make tables
tables = cell(size(phis));
for nChannels_counter = 2 : length(phis)
    values = accuracies_w{nChannels_counter}.accuracies;
    values_rebuilt = accuracies_w{nChannels_counter}.values_rebuilt;
    set_ids = zeros(size(values));
    fly_ids = zeros(size(values));
    
    for set_id = 1 : size(values, 1)
        set_ids(set_id, :) = set_id;
    end
    for fly = 1 : size(values, 2)
        fly_ids(:, fly) = fly;
    end
    
    tables{nChannels_counter} = table(...
        values(:), values_rebuilt(:), set_ids(:), fly_ids(:),...
        'VariableNames', {'value', 'value_rebuilt', 'set_id', 'fly_id'});
    
    tables{nChannels_counter}.fly_id = nominal(tables{nChannels_counter}.fly_id);
    tables{nChannels_counter}.set_id = nominal(tables{nChannels_counter}.set_id);
end

nChannels_counter = 3;

% Build model
model_spec = 'value ~ value_rebuilt + (1|fly_id) + (1|fly_id:set_id)';
model_full = fitlme(tables{nChannels_counter}, model_spec);

model_null_spec = 'value ~ 1 + (1|fly_id) + (1|fly_id:set_id)';
model_null = fitlme(tables{nChannels_counter}, model_null_spec);

compare(model_null, model_full)

%% Regress 3,4ch classification onto 2ch classification (across)

% Make tables
tables = cell(size(phis));
for nChannels_counter = 2 : length(phis)
    values = accuracies_a{nChannels_counter}.accuracies;
    values_rebuilt = accuracies_a{nChannels_counter}.values_rebuilt;
    set_ids = zeros(size(values));
    
    for set_id = 1 : size(values, 1)
        set_ids(set_id, :) = set_id;
    end
    
    tables{nChannels_counter} = table(...
        values(:), values_rebuilt(:), set_ids(:),...
        'VariableNames', {'value', 'value_rebuilt', 'set_id'});
    
    tables{nChannels_counter}.set_id = nominal(tables{nChannels_counter}.set_id);
end

nChannels_counter = 3;

% Build model
model_spec = 'value ~ value_rebuilt + (1|set_id)';
model_full = fitlme(tables{nChannels_counter}, model_spec);

model_null_spec = 'value ~ 1 + (1|set_id)';
model_null = fitlme(tables{nChannels_counter}, model_null_spec);

compare(model_null, model_full)

%% Print figure

% figure_name = 'fig5'; % 'fig5' or 'figS5'
% 
% print(figure_name, '-dsvg'); % SVG
% print(figure_name, '-dpdf', '-bestfit'); % PDF
% print(figure_name, '-dpng'); % PNG

%% Magnitude - Accuracy correlation

nChannels_colours = 'kgb';

flies = (1:13);
condition = 1; % Let's look at correlation between wake value and accuracy

figure;
for nChannels_counter = 1 : length(phis)
    nChannels = nChannels_counter + 1;
    
    % 2.25s trial phi (averaged) vs across-trial classification
    values = mean(mean(phis{nChannels_counter}.phis(:, :, flies, condition, tau), 2), 3);
    accuracy_values = mean(accuracies_w{nChannels_counter}.accuracies(:, flies), 2);
    subplot(1, 2, 1);
    scatter(values, accuracy_values, [nChannels_colours(nChannels_counter) '.']); hold on;
    [r, p] = corr(values, accuracy_values);
    disp(['2.25s/within-trial, ' num2str(nChannels) ' channels: r=' num2str(r) ' p=' num2str(p)]);
    
    % 18s trial phi vs across-fly classification
    values = mean(mean(phis{nChannels_counter}.phis(:, :, flies, condition, tau), 2), 3); % Averaging across trials does nothing (dimension size = 1)
    accuracy_values = mean(accuracies_a{nChannels_counter}.accuracies, 2); % Averaging across flies does nothing (dimension size = 1)
    subplot(1, 2, 2);
    scatter(values, accuracy_values, [nChannels_colours(nChannels_counter) '.']); hold on;
    [r, p] = corr(values, accuracy_values);
    disp(['18s/across-trial, ' num2str(nChannels) ' channels: r=' num2str(r) ' p=' num2str(p)]);
    
end

%% Wake-Anest correlation

nChannels_colours = 'kgb';

flies = (1:13);

figure;
for nChannels_counter = 1 : length(phis)
    nChannels = nChannels_counter + 1;
    
    % 2.25s trial phi (averaged)
    wake = mean(mean(phis{nChannels_counter}.phis(:, :, flies, 1, tau), 2), 3);
    anest = mean(mean(phis{nChannels_counter}.phis(:, :, flies, 2, tau), 2), 3);
    subplot(1, 2, 1);
    scatter(wake, anest, [nChannels_colours(nChannels_counter) '.']); hold on;
    [r, p] = corr(wake, anest);
    disp(['2.25s, ' num2str(nChannels) ' channels: r=' num2str(r) ' p=' num2str(p)]);
    
    % 18s trial phi
    wake = mean(mean(phis_a{nChannels_counter}.phis(:, :, flies, 1, tau), 2), 3); % Averaging across trials does nothing (dimension size = 1)
    anest = mean(mean(phis_a{nChannels_counter}.phis(:, :, flies, 2, tau), 2), 3); % Averaging across trials does nothing (dimension size = 1)
    subplot(1, 2, 2);
    scatter(wake, anest, [nChannels_colours(nChannels_counter) '.']); hold on;
    [r, p] = corr(wake, anest);
    disp(['18s/across-trial, ' num2str(nChannels) ' channels: r=' num2str(r) ' p=' num2str(p)]);
    
end

%% Correlation between 2.25s (averaged) and 18s

nChannels_colours = 'kgb';

flies = (1:13);

condition = 2;

figure;
for nChannels_counter = 1 : length(phis)
    nChannels = nChannels_counter + 1;
    
    % 2.25s trial phi (averaged)
    small = mean(mean(phis{nChannels_counter}.phis(:, :, flies, condition, tau), 2), 3);
    big = mean(mean(phis_a{nChannels_counter}.phis(:, :, flies, condition, tau), 2), 3); % Averaging across trials does nothing (dimension size = 1)
    scatter(small, big, [nChannels_colours(nChannels_counter) '.']); hold on;
    [r, p] = corr(small, big);
    disp(['2.25s/18s, ' num2str(nChannels) ' channels: r=' num2str(r) ' p=' num2str(p)]);
    
end
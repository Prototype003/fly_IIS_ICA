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

%% Load phi values (ACROSS)

[phis_a, measure_string] = phi_load(measure, 1, bin_location);

%% Load accuracy results for phi

results_directory = '../workspace_results/';

if strcmp(measure, 'phi_three')
    % Load ACROSS
    results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
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
    results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
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
    phis_a{nChannels_counter}.channel_sums = zeros(tmp_size);
    phis_a{nChannels_counter}.channel_means = zeros(tmp_size);
    phis_a{nChannels_counter}.set_counters = zeros(tmp_size);
    
    % Accuracies
    tmp_size = size(accuracies_w{nChannels_counter}.accuracies);
    tmp_size(1) = length(channels);
    accuracies_w{nChannels_counter}.channel_sums = zeros(tmp_size);
    accuracies_w{nChannels_counter}.channel_means = zeros(tmp_size);
    accuracies_w{nChannels_counter}.set_counters = zeros(tmp_size);
    accuracies_a{nChannels_counter}.channel_sums = zeros(tmp_size);
    accuracies_a{nChannels_counter}.channel_means = zeros(tmp_size);
    accuracies_a{nChannels_counter}.set_counters = zeros(tmp_size);
    
    % Sum values for each channel (sum across networks which contain the channel)
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                
                phis{nChannels_counter}.channel_sums(channel, :, :, :) = phis{nChannels_counter}.channel_sums(channel, :, :, :) + phis{nChannels_counter}.phis(channel_set, :, :, :, tau);
                phis{nChannels_counter}.set_counters(channel, :, :, :) = phis{nChannels_counter}.set_counters(channel, :, :, :) + 1;
                
                phis_a{nChannels_counter}.channel_sums(channel, :, :, :) = phis_a{nChannels_counter}.channel_sums(channel, :, :, :) + phis_a{nChannels_counter}.phis(channel_set, :, :, :, tau);
                phis_a{nChannels_counter}.set_counters(channel, :, :, :) = phis_a{nChannels_counter}.set_counters(channel, :, :, :) + 1;
                
                accuracies_w{nChannels_counter}.channel_sums(channel, :) = accuracies_w{nChannels_counter}.channel_sums(channel, :) + accuracies_w{nChannels_counter}.accuracies(channel_set, :);
                accuracies_w{nChannels_counter}.set_counters(channel, :) = accuracies_w{nChannels_counter}.set_counters(channel, :) + 1;
                
                accuracies_a{nChannels_counter}.channel_sums(channel, :) = accuracies_a{nChannels_counter}.channel_sums(channel, :) + accuracies_a{nChannels_counter}.accuracies(channel_set, :);
                accuracies_a{nChannels_counter}.set_counters(channel, :) = accuracies_a{nChannels_counter}.set_counters(channel, :) + 1;
                
            end
        end
    end
    
    % Average values
    phis{nChannels_counter}.channel_means = phis{nChannels_counter}.channel_sums ./ phis{nChannels_counter}.set_counters;
    phis_a{nChannels_counter}.channel_means = phis_a{nChannels_counter}.channel_sums ./ phis_a{nChannels_counter}.set_counters;
    accuracies_w{nChannels_counter}.channel_means = accuracies_w{nChannels_counter}.channel_sums ./ accuracies_w{nChannels_counter}.set_counters;
    accuracies_a{nChannels_counter}.channel_means = accuracies_a{nChannels_counter}.channel_sums ./ accuracies_a{nChannels_counter}.set_counters;
    
end

%% Make figure

figure;
set(gcf, 'Position', [0 0 2100/2 600]);
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
% set(gcf, 'RendererMode', 'manual');
% set(gcf, 'Renderer', 'painters');

titles = {'3ch', '4ch'};
value_titles = {'2.25s TPM', '18s TPM'};
class_titles = {'within', 'across'};

ySpacing = 0.2;
xSpacing = 0.05;
xSpacing_wa = 0.1; % spacing between within/across columns

yPortion = 0.8;
xPortion = 0.8;

textbox_width = 0.03;
text_labels = 'abcdefgh';

heights = [1/2 1/2] * yPortion;
widths = ([1 1 1 1] * xPortion ./ 5);

yStarts = (1-yPortion)/2 + fliplr([0 cumsum(fliplr(heights(2:end)))]);
xStarts = (1-xPortion)/2 + (0 : widths(1) : xPortion);
xStarts(2:end) = xStarts(2:end) + xSpacing;
xStarts(3:end) = xStarts(3:end) + xSpacing_wa;
xStarts(4) = xStarts(4) + xSpacing;

measures = {'values', 'class'};
tpm_types = {'nonGlobal', 'global'};

xCounter = 1;
yCounter = 1;
subplot_counter = 1;

colours = 'rb';
shapes = 'ox';

xlims = [...
    0.001 0.004
    0.001 0.004
    0.0005 0.0035
    0.0005 0.0035;
    56.5 60.5;
    56.5 60.5;
    54 63;
    54 63];

for measure_counter = 1 : length(measures)
    measure_type = measures(measure_counter);
    
    xCounter = 1;
    for tpm_type_counter = 1 : length(tpm_types)
        tpm_type = tpm_types{tpm_type_counter};
        
        for nChannels_counter = 2 : length(phis)
            
            subplot(length(measures), length(tpm_types)*2, subplot_counter);
            set(gca, 'Position', [xStarts(xCounter) yStarts(measure_counter)+ySpacing  - (1-yPortion)/3, widths(xCounter), heights(measure_counter)-ySpacing]);
            
            channel_sets = phis{nChannels_counter}.channel_sets;
            if strcmp(measure_type, 'values')
                
                for condition = 1 : 2
                    
                    if strcmp(tpm_type, 'global')
                        values = permute(mean(mean(phis_a{nChannels_counter}.phis(:, :, :, condition, tau), 2), 3), [1 4 2 3]);
                        channel_means = permute(mean(mean(phis_a{1}.channel_means(:, :, :, condition), 2), 3), [1 4 2 3]);
                        
%                         % Can also plot for just one fly here
%                         flies = (1:13);
%                         values = permute(mean(phis_a{nChannels_counter}.phis(:, :, flies, condition, tau), 2), [1 4 2 3]);
%                         channel_means = permute(mean(phis{1}.channel_means(:, :, flies, condition), 2), [1 4 2 3]);
                    else % strcmp(tpm_type, 'nonGlobal')
                        values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, :, condition, tau), 2), 3), [1 4 2 3]);
                        channel_means = permute(mean(mean(phis{1}.channel_means(:, :, :, condition), 2), 3), [1 4 2 3]);
                        
%                         % Can also plot for just one fly here
%                         flies = (1:13);
%                         values = permute(mean(phis{nChannels_counter}.phis(:, :, flies, condition, tau), 2), [1 4 2 3]);
%                         channel_means = permute(mean(phis{1}.channel_means(:, :, flies, condition), 2), [1 4 2 3]);
                    end
                    
                    % Build values from average of single channel average values
                    values_rebuilt = zeros(size(values));
                    for channel_set = 1 : size(channel_sets, 1)
                        channels = channel_sets(channel_set, :);
                        values_rebuilt(channel_set) = mean(channel_means(channels));
                    end
                    
                    % Plot
                    scatter(values_rebuilt, values, [], colours(condition), shapes(condition)); hold on;
                    
                    % Slope
                    [bs] = [ones(size(values_rebuilt)) values_rebuilt] \ values;
                    disp([measure_type tpm_type num2str(nChannels_counter) 'c' num2str(condition) ' slope:' num2str(bs(2)) ' intercept:' num2str(bs(1))]);
                    
                    % Correlation
                    [r, p] = corr(values_rebuilt, values);
                    disp([measure_type tpm_type num2str(nChannels_counter) 'c' num2str(condition) ' r:' num2str(r) ' p:' num2str(p)]);
                end
                
                title([num2str(nChannels_counter+1) 'ch']);
                
                if nChannels_counter == 2
                    ylabel(measure_string);
                    xlabel(['mean 2ch ' measure_string]);
                end
                
                handle = gca;
                handle.XRuler.Exponent = 0;
                
            else %strcmp(measure_type, 'class')
                
                if strcmp(tpm_type, 'global')
                    values = mean(accuracies_a{nChannels_counter}.accuracies, 2);
                    channel_means = mean(accuracies_a{1}.channel_means, 2);
                else % strcmp(tpm_type, 'nonGlobal')
                    values = mean(accuracies_w{nChannels_counter}.accuracies, 2);
                    channel_means = mean(accuracies_w{1}.channel_means, 2);
                end
                
                % Build values from average of single channel average values
                values_rebuilt = zeros(size(values));
                for channel_set = 1 : size(channel_sets, 1)
                    channels = channel_sets(channel_set, :);
                    values_rebuilt(channel_set) = mean(channel_means(channels));
                end
                
                % Plot
                scatter(values_rebuilt, values, [], 'k.'); hold on;
                
                % Correlation
                [r, p] = corr(values_rebuilt, values);
                disp([measure_type tpm_type num2str(nChannels_counter) ' r:' num2str(r) ' p:' num2str(p)]);
                
                if nChannels_counter == 2
                    ylabel('class. acc. %');
                    xlabel('mean 2ch acc. %');
                end
                
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
            
            % Add column group heading
            if nChannels_counter == 2
                % Phi values
                axes('Visible', 'off', 'Position', [(xStarts(xCounter)+widths(xCounter)+xStarts(xCounter+1))/2 yStarts(1)+heights(1), 0.3 0.05]); axis_defaults(gca);
                text(0, 0, value_titles{tpm_type_counter}, 'FontSize', fontsize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                % Classification
                axes('Visible', 'off', 'Position', [(xStarts(xCounter)+widths(xCounter)+xStarts(xCounter+1))/2 yStarts(2)+heights(2), 0.3 0.05]); axis_defaults(gca);
                text(0, 0, class_titles{tpm_type_counter}, 'FontSize', fontsize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            end
            
            xCounter = xCounter + 1;
            subplot_counter = subplot_counter + 1;
        end
        
    end
    
    yCounter = yCounter + 1;
end

%% Print figure

% figure_name = 'fig5';
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
    values = mean(mean(phis_a{nChannels_counter}.phis(:, :, flies, condition, tau), 2), 3); % Averaging across trials does nothing (dimension size = 1)
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
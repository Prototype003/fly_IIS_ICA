%% DESCRIPTION

%{

Figure 4 - classification boxplots / bar + center/distance plots (classification)
1 row + 3 x 2 (2/3/4ch by within/across)
OR
1 column (horizontal box/plot) + 3 x 2 (2/3/4ch by within/across

Run fig4_tscore_values.m to get values to plot

%}

%% Setup

load('fig4_tscore.mat');
measure_string = '\Phi^{3.0}';

%% Setup

measure = 'phi_three';
tau = 1; % 1 = 4ms; 2 = 8ms; 3 = 16ms
if tau == 1
    tau_string = '4';
elseif tau == 2
    tau_string = '8';
elseif tau == 3
    tau_string = '16';
end

freq_range_w = (1:42); %(1:83); % corresponding to ~5Hz and ~10Hz, check the 'frequencies' vector
freq_range_a = (1:329); %(1:329)=0-5Hz; There are more frequency bins for the single large trial
freq_range_string = '0-5Hz'; %'0-10Hz';

fontsize = 11; % Used for drawing label letters

bin_location = '../';
addpath(bin_location);

results_directory = [bin_location 'workspace_results/'];

% %% Data for boxplots (within)
% 
% % measure_accuracies - this held classification accuracy, but now will hold t-scores
% 
% % Power
% results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification.mat';
% load([results_directory results_filename]);
% powers = permute(mean(powers(freq_range_w, :, :, :, :), 1), [2 3 5 4 1]); % mean across frequency range; trials x channels x flies x conditions
% scores = zeros(size(powers, 3), size(powers, 2)); % flies x channels
% for channel = 1 : size(powers, 2)
%     for fly = 1 : size(powers, 3)
%         % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
%         [h, p, ci, stats] = ttest2(powers(:, channel, fly, 1), powers(:, channel, fly, 2));
%         %[p, h, stats] = ranksum(powers(:, channel, fly, 1), powers(:, channel, fly, 2));
%         scores(fly, channel) = stats.tstat; % stats.ranksum; stats.tstat;
%     end
% end
% % Average scores across channels
% scores_mean = mean(scores, 2);
% % Add to boxplot data structure
% measure_accuracies_w = scores_mean;
% measure_groups_w = zeros(size(measure_accuracies_w)) + 1;
% 
% 
% 
% % Coherence
% results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';
% load([results_directory results_filename]);
% coherencies = permute(mean(coherencies(freq_range_w, :, :, :, :), 1), [2 3 5 4 1]); % mean across frequency range; trials x ch-pairs x flies x conditions
% scores = zeros(size(coherencies, 3), size(coherencies, 2)); % flies x ch-pairs
% for pair = 1 : size(coherencies, 2)
%     for fly = 1 : size(coherencies, 3)
%         % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
%         [h, p, ci, stats] = ttest2(coherencies(:, pair, fly, 1), coherencies(:, pair, fly, 2));
%         %[p, h, stats] = ranksum(coherencies(:, pair, fly, 1), coherencies(:, pair, fly, 2));
%         scores(fly, pair) = stats.tstat; % stats.ranksum; stats.tstat;
%     end
% end
% % Average scores across pairs
% scores_mean = mean(scores, 2);
% % Add to boxplot data structure
% measure_accuracies_w = [measure_accuracies_w; scores_mean];
% measure_groups_w = [measure_groups_w; zeros(size(scores_mean)) + 2];
% 
% 
% 
% % Phi-three
% [phis, measure_string] = phi_load('phi_three', 0, bin_location);
% group_counter = 3;
% for nChannels_counter = 1 : length(phis)
%     values = phis{nChannels_counter}.phis(:, :, :, :, tau);
%     phis{nChannels_counter}.tscores = zeros(size(values, 3), size(values, 1)); % flies x ch-sets
%     for set_counter = 1 : size(values, 1)
%         for fly = 1 : size(values, 3)
%             % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
%             [h, p, ci, stats] = ttest2(values(set_counter, :, fly, 1), values(set_counter, :, fly, 2));
%             %[p, h, stats] = ranksum(values(set_counter, :, fly, 1), values(set_counter, :, fly, 2));
%             phis{nChannels_counter}.tscores(fly, set_counter) = stats.tstat; %stats.ranksum; stats.tstat;
%         end
%     end
%     % Average scores across sets
%     scores_mean = mean(phis{nChannels_counter}.tscores, 2);
%     % Add to boxplot data structure
%     measure_accuracies_w = [measure_accuracies_w; scores_mean];
%     measure_groups_w = [measure_groups_w; zeros(size(scores_mean)) + group_counter];
%     group_counter = group_counter + 1;
% end
% % Store for other plots (store in accuracies - in order to reuse code)
% if strcmp(measure, 'phi_three')
%     accuracies_w = cell(size(phis));
%     for nChannels_counter = 1 : length(phis)
%         accuracies_w{nChannels_counter}.accuracies = permute(phis{nChannels_counter}.tscores, [2 1]);
%         accuracies_w{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
%     end
% end
% 
% 
% 
% % Phi-star
% [phis, measure_string] = phi_load('phi_star', 0, bin_location);
% for nChannels_counter = 1 : length(phis)
%     values = phis{nChannels_counter}.phis(:, :, :, :, tau);
%     phis{nChannels_counter}.tscores = zeros(size(values, 3), size(values, 1)); % flies x ch-sets
%     for set_counter = 1 : size(values, 1)
%         for fly = 1 : size(values, 3)
%             % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
%             [h, p, ci, stats] = ttest2(values(set_counter, :, fly, 1), values(set_counter, :, fly, 2));
%             %[p, h, stats] = ranksum(values(set_counter, :, fly, 1), values(set_counter, :, fly, 2));
%             phis{nChannels_counter}.tscores(fly, set_counter) = stats.tstat; %stats.ranksum; stats.tstat;
%         end
%     end
%     % Average scores across sets
%     scores_mean = mean(phis{nChannels_counter}.tscores, 2);
%     % Add to boxplot data structure
%     measure_accuracies_w = [measure_accuracies_w; scores_mean];
%     measure_groups_w = [measure_groups_w; zeros(size(scores_mean)) + group_counter];
%     group_counter = group_counter + 1;
% end
% % Store for other plots (store in accuracies - in order to reuse code)
% if strcmp(measure, 'phi_star')
%     accuracies_w = cell(size(phis));
%     for nChannels_counter = 1 : length(phis)
%         accuracies_w{nChannels_counter}.accuracies = permute(phis{nChannels_counter}.tscores, [2 1]);
%         accuracies_w{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
%     end
% end
% 
% %% Data for boxplots (across)
% 
% % Power
% results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across1.mat';
% load([results_directory results_filename]);
% powers = permute(mean(powers(freq_range_a, :, :, :, :), 1), [5 3 4 1 2]); % mean across frequency range; flies x channels x conditions
% scores = zeros(size(powers, 2), 1); % channels
% for channel = 1 : size(powers, 2)
%     % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
%     [h, p, ci, stats] = ttest(powers(:, channel, 1), powers(:, channel, 2));
%     %[p, h, stats] = signrank(powers(:, channel, 1), powers(:, channel, 2));
%     scores(channel) = stats.tstat; % stats.signedrank; stats.tstat;
% end
% % Add to boxplot data structure
% measure_accuracies_a = scores;
% measure_groups_a = zeros(size(measure_accuracies_a)) + 1;
% 
% 
% 
% % Coherence
% results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification_across1.mat';
% load([results_directory results_filename]);
% coherencies = permute(mean(coherencies(freq_range_a, :, :, :, :), 1), [5 3 4 1 2]); % mean across frequency range; flies x ch-pairs x conditions
% scores = zeros(size(coherencies, 2), 1); % ch-pairs
% for pair = 1 : size(coherencies, 2)
%     % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
%     [h, p, ci, stats] = ttest(coherencies(:, pair, 1), coherencies(:, pair, 2));
%     %[p, h, stats] = signrank(coherencies(:, pair, 1), coherencies(:, pair, 2));
%     scores(pair) = stats.tstat; % stats.signedrank; stats.tstat;
% end
% % Add to boxplot data structure
% measure_accuracies_a = [measure_accuracies_a; scores];
% measure_groups_a = [measure_groups_a; zeros(size(scores)) + 2];
% 
% 
% 
% % Phi-three
% [phis, measure_string] = phi_load('phi_three', 1, bin_location);
% group_counter = 3;
% for nChannels_counter = 1 : length(phis)
%     values = permute(phis{nChannels_counter}.phis(:, :, :, :, tau), [3 1 4 2 5]); % flies x sets x conditions
%     phis{nChannels_counter}.tscores = zeros(size(values, 2), 1); % ch-sets
%     for set_counter = 1 : size(values, 2)
%             % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
%             [h, p, ci, stats] = ttest(values(:, set_counter, 1), values(:, set_counter, 2));
%             %[p, h, stats] = signrank(values(:, set_counter, 1), values(:, set_counter, 2));
%             phis{nChannels_counter}.tscores(set_counter) = stats.tstat; %stats.signedrank; stats.tstat;
%     end
%     % Add to boxplot data structure
%     scores = phis{nChannels_counter}.tscores;
%     measure_accuracies_a = [measure_accuracies_a; scores];
%     measure_groups_a = [measure_groups_a; zeros(size(scores)) + group_counter];
%     group_counter = group_counter + 1;
% end
% % Store for other plots (store in accuracies - in order to reuse code)
% if strcmp(measure, 'phi_three')
%     accuracies_a = cell(size(phis));
%     for nChannels_counter = 1 : length(phis)
%         accuracies_a{nChannels_counter}.accuracies = phis{nChannels_counter}.tscores;
%         accuracies_a{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
%     end
% end
% 
% 
% 
% % Phi-star
% [phis, measure_string] = phi_load('phi_star', 1, bin_location);
% for nChannels_counter = 1 : length(phis)
%     values = permute(phis{nChannels_counter}.phis(:, :, :, :, tau), [3 1 4 2 5]); % flies x sets x conditions
%     phis{nChannels_counter}.tscores = zeros(size(values, 2), 1); % ch-sets
%     for set_counter = 1 : size(values, 2)
%             % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
%             [h, p, ci, stats] = ttest(values(:, set_counter, 1), values(:, set_counter, 2));
%             %[p, h, stats] = signrank(values(:, set_counter, 1), values(:, set_counter, 2));
%             phis{nChannels_counter}.tscores(set_counter) = stats.tstat; %stats.signedrank; stats.tstat;
%     end
%     % Add to boxplot data structure
%     scores = phis{nChannels_counter}.tscores;
%     measure_accuracies_a = [measure_accuracies_a; scores];
%     measure_groups_a = [measure_groups_a; zeros(size(scores)) + group_counter];
%     group_counter = group_counter + 1;
% end
% % Store for other plots (store in accuracies - in order to reuse code)
% if strcmp(measure, 'phi_star')
%     accuracies_a = cell(size(phis));
%     for nChannels_counter = 1 : length(phis)
%         accuracies_a{nChannels_counter}.accuracies = phis{nChannels_counter}.tscores;
%         accuracies_a{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
%     end
% end

%% Make figure

figure;
set(gcf, 'Position', [0 0 2100/4 750]);
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
% set(gcf, 'RendererMode', 'manual');
% set(gcf, 'Renderer', 'painters');

colormap('jet');
subplot_counter = 1;
titles = {'within', 'across'};
cbar_titles = {...
    [measure_string '_{W}'],...
    [measure_string '_{A}'],...
    [measure_string '_{W} / ' measure_string '_{A}']...
    };

nChannels_yAspects = [1 2 4] ./ 2;

ySpacing = 0.05;
xSpacing = 0.03;

yPortion = 0.7;
xPortion = 0.7;

textbox_width = 0.03;
text_labels = 'aabcdefg';

cbar_width = 0.01;
cbar_spacing = 0.08;

heights = [7/18 3/18 4/18 7/18] * yPortion; % y-axis doubles each time, x + 2x + 4x = 1
widths = ([1 1] * xPortion ./ 2) - cbar_width;

yStarts = (1-yPortion)/2 + fliplr([0 cumsum(fliplr(heights(2:end)))]);
xStarts = (1-xPortion)/2 + (0 : widths(1) : 1);
xStarts(2) = xStarts(2) + xSpacing;
xStarts(3) = xStarts(3) + cbar_width + xSpacing + cbar_spacing;

% Boxplot
box_widths = 0.2;
condition_offsets = [-0.15 0.15];
box_lineWidths = 1;
box_colors = 'kkbbbkkk';
whisker_length = 1000; % Default 1.5, this determines what points are treated as outliers

% Plot
subplot(length(accuracies_a)+1, 2, [1 2]);
boxplot(measure_accuracies_w, measure_groups_w, 'Positions', measure_groups_w + condition_offsets(1), 'Widths', box_widths, 'MedianStyle', 'target', 'Colors', box_colors, 'Whisker', whisker_length); hold on;
boxplot(measure_accuracies_a, measure_groups_a, 'Positions', measure_groups_a + condition_offsets(2), 'Widths', box_widths, 'MedianStyle', 'target', 'Colors', box_colors, 'Whisker', whisker_length);
% % Means
% for measure = 1 : max(measure_groups_w)
%     scatter(measure + condition_offsets(1), mean(measure_accuracies_w(measure_groups_w==measure)), 'bo');
%     scatter(measure + condition_offsets(2), mean(measure_accuracies_a(measure_groups_a==measure)), 'bo');
% end
%ylim([0 100]);

% Position
set(gca, 'Position', [xStarts(1) yStarts(1)+ySpacing*2  - (1-yPortion)/4, sum(widths)+xSpacing+cbar_width, heights(1)-ySpacing]);

boxes = findobj(gca, 'Tag', 'Box'); % order is presumably from most recent drawn object (boxes(1)) to oldest (boxes(end))

% Set across-flies to dotted line
for box = 1 : length(box_colors)
    boxes(box).LineStyle = ':';
end

% Set line width for all boxes
for box = 1 : length(box_colors)*2
    boxes(box).LineWidth = box_lineWidths;
end

%grid on
grid minor
axis_defaults(gca);
measure_names = {'  P\newline1ch', '  C\newline2ch', ' \Phi^{3.0}\newline 2ch', ' \Phi^{3.0}\newline 3ch', ' \Phi^{3.0}\newline 4ch', ' \Phi*\newline2ch', ' \Phi*\newline3ch', ' \Phi*\newline4ch'};
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'XTick', (1:length(box_colors)));
set(gca, 'XTickLabels', measure_names);
ylabel('t');

% Add letter label
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(subplot_counter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');



% center/distance plots
subplot_counter = 3;
xCounter = 1;
yCounter = 2; % Because the first one was used by the boxplot
for nChannels_counter = 1 : length(accuracies_a)
    nChannels = nChannels_counter + 1;
    channel_sets = accuracies_a{nChannels_counter}.channel_sets;
    
    % Get values to plot
    plot_values = zeros(length(accuracies_a{nChannels_counter}.accuracies), 2);
    plot_values(:, 1) = mean(accuracies_w{nChannels_counter}.accuracies, 2);
    plot_values(:, 2) = accuracies_a{nChannels_counter}.accuracies;
    
    % Set center / set path distance mapping
    centers = mean(channel_sets, 2); % Mean across channels in each set
    distances = channel_set_distances(channel_sets);
    
    % Find minimum delta among map values
    centers_deltas = pdist2(unique(centers), unique(centers));
    centers_deltas(centers_deltas==0) = NaN; % remove 0 distances
    centers_delta = min(centers_deltas(:));
    distances_deltas = pdist2(unique(distances), unique(distances));
    distances_deltas(distances_deltas==0) = NaN; % remove 0 distances
    distances_delta = min(distances_deltas(:));
    
    % Create space mapping (used for indexing into the mapped space)
    %     centers_axis = (min(centers) - (2*centers_delta) : centers_delta : max(centers) + (2*centers_delta)); % - and + delta is for padding
    %     distances_axis = (min(distances) - (2*distances_delta) : distances_delta : max(distances) + (2*distances_delta)); % - and + delta is for padding
    %     [centers_map, distances_map] = meshgrid(centers_axis, distances_axis);
    centers_axis = (min(channel_sets(:)) : centers_delta : max(channel_sets(:))+centers_delta); % range is from first channel to last channel; +delta is for padding (because of how pcolor works)
    distances_axis = (min(distances) : distances_delta : max(distances)+distances_delta); % range is from min to max
    [centers_map, distances_map] = meshgrid(centers_axis, distances_axis);
    
    % Map channel sets to 2D matrix (and average overlaps)
    for value_type = 1 : size(plot_values, 2)
        
        % Create mapped space
        values_map = zeros(size(centers_map)); % Will sum all values with the same coordinates
        values_map_counter = zeros(size(centers_map)); % Keeps count in each coordinate as to how many values have that coordinate
        
        % Populate mapped space
        for value_counter = 1 : size(plot_values, 1)
            x = find(abs(centers_axis - centers(value_counter)) < 0.00001, 1); % This gives the mapped x location
            y = find(distances_axis == distances(value_counter), 1); % This gives the mapped y location
            values_map(y, x) = values_map(y, x) + plot_values(value_counter, value_type); % matrix is (rows, columns), corresponding to (y, x)
            values_map_counter(y, x) = values_map_counter(y, x) + 1;
        end
        
        values_map_plot = values_map ./ values_map_counter; % We will use the NaN values for the black background (and maybe interpolation)
        
        % Interpolate values (linearly - each NaN will turn into the average of the cells immediately adjacent to it, excluding diagonals)
        min_value = min(values_map_plot(:)); % We will remove interpolations which are less than the original min value in the plot (to avoid 'ghost/shadow interpolations');
        values_adjacent = zeros(4, 1);
        for y = 2 : size(values_map_plot, 1) - 1
            for x = 2 : size(values_map_plot, 2) - 1
                if isnan(values_map_plot(y, x))
                    values_adjacent(1) = values_map_plot(y, x-1);
                    values_adjacent(2) = values_map_plot(y, x+1);
                    values_adjacent(3) = values_map_plot(y-1, x);
                    values_adjacent(4) = values_map_plot(y+1, x);
                    if sum(~isnan(values_adjacent)) > 2
                        %values_adjacent(isnan(values_adjacent)) = 0; % Deal with NaNs
                        %values_adjacent(values_adjacent < min_value) = 0; % Don't include 'shadow' interpolations
                        interpolated = mean(values_adjacent(~isnan(values_adjacent)));
                        if interpolated == 0
                            values_map_plot(y, x) = NaN; % Place NaN instead of 0 to avoid affecting the color scaling
                        else
                            values_map_plot(y, x) = interpolated;
                        end
                    end
                end
            end
        end
        
        if value_type == 1
            % Create storage matrix for all value types
            values_map_plots = zeros([size(values_map_plot) size(plot_values, 2)]);
        end
        values_map_plots(:, :, value_type) = values_map_plot;
    end
    
    % Consistent colorbar within set size
    clim = round([min(values_map_plots(:)) max(values_map_plots(:))]);
    %clim = [42 80];
    
    % Plot
    xCounter = 1;
    for value_type = 1 : size(plot_values, 2)
        subplot(length(accuracies_a)+1, size(plot_values, 2), subplot_counter);
        set(gca, 'Position', [xStarts(xCounter) yStarts(nChannels_counter+1)+ySpacing  - (1-yPortion)/3, widths(xCounter), heights(nChannels_counter+1)-ySpacing]);
        
        plot = pcolor(values_map_plots(:, :, value_type));
        
        axis_defaults(gca);
        
        set(gca, 'TickDir', 'out');
        set(gca, 'Box', 'on');
        
        if value_type < 3
            
            % Set shared colour axis
            caxis(clim);
            
        end
        if value_type > 1
            
            % Add colorbar (and reposition)
            ax_pos = get(gca, 'Position');
            c = colorbar;
            set(c, 'Position', [ax_pos(1)+ax_pos(3) ax_pos(2) cbar_width ax_pos(4)]); % Reposition colourbar directly next to plot
            set(c, 'Box', 'on');
            c.Ruler.Exponent = 0; % Get rid of scientific notation
            %set(c, 'Ticks', [1, ], 'TickLabels', [])% label min, mid, max values
            clims = get(c, 'Limits');
            set(c, 'Ticks', linspace(clims(1), clims(2), 3)); % Set ticks to be min, mid, and max
%             if nChannels_counter == 1 % Set decimal places to display
%                 c.Ruler.TickLabelFormat = '%.4f';
%             else
%                 c.Ruler.TickLabelFormat = '%.3f';
%             end
%             if value_type == 3
%                 c.Ruler.TickLabelFormat = '%.2f';
%             end
            set(gca, 'Position', ax_pos); % Reposition plot to original specs
            
            % Hide YTicks
            %set(gca, 'YTick', []);
            set(gca, 'YTickLabel', []);
        end
        
        if nChannels_counter == 1
            title(titles{value_type});
            if value_type == 2
                set(get(c, 'Title'), 'String', '%')
            end
        end
        
        set(gca, 'XTick', [1.5 ceil(length(centers_axis)/2)+0.5 length(centers_axis)-0.5]);
        set(gca, 'YTick', [1.5 ceil(length(distances_axis)/2)+0.5 length(distances_axis)-0.5]);
        
        if value_type == 1
            set(gca, 'YTickLabel', [min(distances) median(distances) max(distances)]);
        end
        
        if nChannels_counter == 3 && value_type == 1
            set(gca, 'XTickLabel', [1 median(centers) 15]);
        else
            set(gca, 'XTickLabel', []);
        end
        
        % X and Y axis labels
        if nChannels_counter == 3 && value_type == 2
            xlabel('set centre');
        end
        if nChannels_counter == 2 && value_type == 1
            ylabel('set path distance');
        end
        
        set(gca, 'color', [0 0 0]); % black background
        set(plot, 'EdgeColor', 'none'); % remove grid outline
        
        % Add letter label to plot
        ax_pos = get(gca, 'Position');
        axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
        text(0, 0.2, text_labels(subplot_counter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');
        
        subplot_counter = subplot_counter + 1;
        xCounter = xCounter + 1;
    end
end

%% Print figure

% figure_name = 'fig4';
% 
% print(figure_name, '-dsvg'); % SVG
% print(figure_name, '-dpdf', '-bestfit'); % PDF
% print(figure_name, '-dpng'); % PNG
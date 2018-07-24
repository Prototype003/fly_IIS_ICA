%% DESCRIPTION

%{

Figure 4 - classification boxplots / bar + center/distance plots (classification)
1 row + 3 x 2 (2/3/4ch by within/across)
OR
1 column (horizontal box/plot) + 3 x 2 (2/3/4ch by within/across


%}

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

% 0 for tpms/covariances constructed per trial, and for within fly analyses
% 1 for tpms/covariances constructed across trials, and for across fly analyses
global_tpm = 0;

fontsize = 11; % Used for drawing label letters

bin_location = '../';
addpath(bin_location);

%% Load files
% 
% [phis, measure_string] = phi_load(measure, global_tpm, bin_location);

%% Get t-scores
% This is for within flies (t-score across trials)

for nChannels_counter = 1 : length(phis)
    mat_size = size(phis{nChannels_counter}.phis(:, :, :, :, tau));
    mat_size = mat_size([1 3]);
    phis{nChannels_counter}.tscores = zeros(mat_size);
    for set_counter = 1 : mat_size(1)
        for fly = 1 : mat_size(2)
            % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
            %[h, p, ci, stats] = ttest2(phis{nChannels_counter}.phis(set_counter, :, fly, 1, tau), phis{nChannels_counter}.phis(set_counter, :, fly, 2, tau));
            [p, h, stats] = ranksum(phis{nChannels_counter}.phis(set_counter, :, fly, 1, tau), phis{nChannels_counter}.phis(set_counter, :, fly, 2, tau));
            phis{nChannels_counter}.tscores(set_counter, fly) = stats.ranksum; %stats.signrank; stats.tscore;
            
            % Cohen's D
            %phis{nChannels_counter}.tscores(set_counter, fly) =...
             %   (mean(phis{nChannels_counter}.phis(set_counter, :, fly, 1, tau)) - mean(phis{nChannels_counter}.phis(set_counter, :, fly, 2, tau)))...
              %  / sqrt((std(phis{nChannels_counter}.phis(set_counter, :, fly, 1, tau))+std(phis{nChannels_counter}.phis(set_counter, :, fly, 2, tau)))/2);
        end
    end
end

%% Center-distance plot
% This is for within flies - plots average t-scores

figure;
colormap('jet');
subplot_counter = 1;
for nChannels_counter = 1 : length(phis)
    nChannels = phis{nChannels_counter}.nChannels;
    channel_sets = phis{nChannels_counter}.channel_sets;
    channel_sets = double(channel_sets);
    
    % Get values to plot
    plot_values = mean(phis{nChannels_counter}.tscores, 2);
    
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
    
    % Create mapped space
    values_map = zeros(size(centers_map)); % Will sum all values with the same coordinates
    values_map_counter = zeros(size(centers_map)); % Keeps count in each coordinate as to how many values have that coordinate
    
    % Map channel sets to 2D matrix (and average overlaps)
    for value_type = 1 : size(plot_values, 2)
        
        % Populate mapped space
        for value_counter = 1 : size(plot_values, 1)
            x = find(abs(centers_axis - centers(value_counter)) < 0.00001, 1); % This gives the mapped x location
            y = find(distances_axis == distances(value_counter), 1); % This gives the mapped y location
            values_map(y, x) = values_map(y, x) + plot_values(value_counter, value_type); % matrix is (rows, columns), corresponding to (y, x)
            values_map_counter(y, x) = values_map_counter(y, x) + 1;
        end
        
        values_map_plot = values_map ./ values_map_counter; % We will use the NaN values for the black background (and maybe interpolation)
        
%         % Interpolate values (linearly - each NaN will turn into the average of the cells immediately adjacent to it, excluding diagonals)
%         min_value = min(values_map_plot(:)); % We will remove interpolations which are less than the original min value in the plot (to avoid 'ghost/shadow interpolations');
%         values_adjacent = zeros(4, 1);
%         for y = 2 : size(values_map_plot, 1) - 1
%             for x = 2 : size(values_map_plot, 2) - 1
%                 if isnan(values_map_plot(y, x))
%                     values_adjacent(1) = values_map_plot(y, x-1);
%                     values_adjacent(2) = values_map_plot(y, x+1);
%                     values_adjacent(3) = values_map_plot(y-1, x);
%                     values_adjacent(4) = values_map_plot(y+1, x);
%                     if sum(~isnan(values_adjacent)) > 2
%                         %values_adjacent(isnan(values_adjacent)) = 0; % Deal with NaNs
%                         %values_adjacent(values_adjacent < min_value) = 0; % Don't include 'shadow' interpolations
%                         interpolated = mean(values_adjacent(~isnan(values_adjacent)));
%                         if interpolated == 0
%                             values_map_plot(y, x) = NaN; % Place NaN instead of 0 to avoid affecting the color scaling
%                         else
%                             values_map_plot(y, x) = interpolated;
%                         end
%                     end
%                 end
%             end
%         end
        
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
        subplot(length(phis), 1, subplot_counter);
        %set(gca, 'Position', [xStarts(xCounter) yStarts(nChannels_counter+1)+ySpacing  - (1-yPortion)/3, widths(xCounter), heights(nChannels_counter+1)-ySpacing]);
        
        plot = pcolor(values_map_plots(:, :, value_type)); colorbar;
        
%         axis_defaults(gca);
%         
%         set(gca, 'TickDir', 'out');
%         set(gca, 'Box', 'on');
%         
%         if value_type < 3
%             
%             % Set shared colour axis
%             caxis(clim);
%             
%         end
%         if value_type > 1
%             
%             % Add colorbar (and reposition)
%             ax_pos = get(gca, 'Position');
%             c = colorbar;
%             set(c, 'Position', [ax_pos(1)+ax_pos(3) ax_pos(2) cbar_width ax_pos(4)]); % Reposition colourbar directly next to plot
%             set(c, 'Box', 'on');
%             c.Ruler.Exponent = 0; % Get rid of scientific notation
%             %set(c, 'Ticks', [1, ], 'TickLabels', [])% label min, mid, max values
%             clims = get(c, 'Limits');
%             set(c, 'Ticks', linspace(clims(1), clims(2), 3)); % Set ticks to be min, mid, and max
% %             if nChannels_counter == 1 % Set decimal places to display
% %                 c.Ruler.TickLabelFormat = '%.4f';
% %             else
% %                 c.Ruler.TickLabelFormat = '%.3f';
% %             end
% %             if value_type == 3
% %                 c.Ruler.TickLabelFormat = '%.2f';
% %             end
%             set(gca, 'Position', ax_pos); % Reposition plot to original specs
%             
%             % Hide YTicks
%             %set(gca, 'YTick', []);
%             set(gca, 'YTickLabel', []);
%         end
%         
%         if nChannels_counter == 1
%             title(titles{value_type});
%             if value_type == 2
%                 set(get(c, 'Title'), 'String', '%')
%             end
%         end
%         
%         set(gca, 'XTick', [1.5 ceil(length(centers_axis)/2)+0.5 length(centers_axis)-0.5]);
%         set(gca, 'YTick', [1.5 ceil(length(distances_axis)/2)+0.5 length(distances_axis)-0.5]);
%         
%         if value_type == 1
%             set(gca, 'YTickLabel', [min(distances) median(distances) max(distances)]);
%         end
%         
%         if nChannels_counter == 3 && value_type == 1
%             set(gca, 'XTickLabel', [1 median(centers) 15]);
%         else
%             set(gca, 'XTickLabel', []);
%         end
%         
%         % X and Y axis labels
%         if nChannels_counter == 3 && value_type == 2
%             xlabel('set centre');
%         end
%         if nChannels_counter == 2 && value_type == 1
%             ylabel('set path distance');
%         end
%         
        set(gca, 'color', [0 0 0]); % black background
        set(plot, 'EdgeColor', 'none'); % remove grid outline
%         
%         % Add letter label to plot
%         ax_pos = get(gca, 'Position');
%         axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
%         text(0, 0.2, text_labels(subplot_counter), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');
%         
        subplot_counter = subplot_counter + 1;
%         xCounter = xCounter + 1;
    end
end

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

% center/distance plots
subplot_counter = 3;
xCounter = 1;
yCounter = 2; % Because the first one was used by the boxplot
for nChannels_counter = 1 : length(phis)
    nChannels = phis{nChannels_counter}.nChannels;
    channel_sets = phis{nChannels_counter}.channel_sets;
    
    % Get values to plot
    values = mean(accuracies{nChannels_counter}.accuracies, 2);
    plot_values = zeros(size(values, 1), 2);
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
    
    % Create mapped space
    values_map = zeros(size(centers_map)); % Will sum all values with the same coordinates
    values_map_counter = zeros(size(centers_map)); % Keeps count in each coordinate as to how many values have that coordinate
    
    % Map channel sets to 2D matrix (and average overlaps)
    for value_type = 1 : size(plot_values, 2)
        
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
        subplot(length(phis)+1, size(plot_values, 2), subplot_counter);
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

figure_name = 'fig4';

print(figure_name, '-dsvg'); % SVG
print(figure_name, '-dpdf', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG
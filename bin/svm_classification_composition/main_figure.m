%% Description

%{

Make summary figure

Average across sets or across flies?

Correlation between measures - does phi perform better when power also
performs better?

Include results from 15ch classification

center-distance map - is it still worth including?

%}

%% Setup

measures = {'4ch_medianSplit_power', '4ch_medianSplit_coherence', '2ch_phi3Composition_unpart', '3ch_phi3Composition_unpart', '4ch_phi3Composition_unpart', '4ch_phiStarComposition', '4ch_phiStarGaussianComposition'};
measure_labels = {'P', 'C', '2ch\phi3comp', '3ch\phi3comp', '4ch\phi3comp', '\Phi*comp', '\Phi*gcomp'};

clim = [0 1];

%% Load within-fly results (composition classification)

results_location = 'results/';

results = cell(size(measures));
results_mat = []; % set x fly x measure

class_type = 'within';

for measure = 1 : length(measures)
    filename = [measures{measure} '_svm_' class_type '.mat'];
    results{measure} = load([results_location filename]);
    results{measure}.accuracy = squeeze(results{measure}.accuracy); % Power/coherence have a leading singleton dimension
    results_mat = cat(3, results_mat, results{measure}.accuracy);
end

% Average across flies
results_comp = permute(mean(results_mat, 2), [1 3 2]);

%% Load across-fly results (composition classification)

results_location = 'results/';

results = cell(size(measures));
results_mat = []; % set x measure

class_type = 'across';

for measure = 1 : length(measures)
    filename = [measures{measure} '_svm_' class_type '.mat'];
    results{measure} = load([results_location filename]);
    results{measure}.accuracy = permute(results{measure}.accuracy, [2 1]);
    results_mat = cat(2, results_mat, results{measure}.accuracy);
end

% Join with other results
results_comp = cat(3, results_comp, results_mat);

%% Load within-fly results (15ch classification)

results_location = '../svm_classification_15ch/results/';
measures = {'medianSplit_power', 'medianSplit_coherence', 'phi3', 'phi3CompositionOverlapping', 'PhiStar', 'PhiStarComposition', 'PhiStarGaussian', 'PhiStarGaussianComposition'};
measure_labels = {'P', 'C', '\Phi^{3}', '\phicomp', '\Phi*', '\Phi*comp', '\Phi*g', '\Phi*gcomp'};

results = cell(1, 2);

class_type = 'within';

for measure = 1 : length(measures)
    results{measure} = load([results_location measures{measure} '_svm_' class_type '.mat']);
end

% Average accuracies across flies
within = [];
for measure = 1 : length(results)
    within = cat(2, within, results{measure}.accuracy);
end

results_overlap = within;

%% Load across-fly results (15ch classification)

%% Channel-set centre-distance mapping
% TODO: turn into a function

plot_measure = 3;

for plot_measure = 1 : length(measures)

addpath('../');

nChannels = 4;
channel_sets = nchoosek((1:15), 4);

% Get values to plot
plot_values = zeros(size(results_comp, 1), 2);
plot_values(:, 1) = results_comp(:, plot_measure, 1);
plot_values(:, 2) = results_comp(:, plot_measure, 2);

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

% Make figure

figure;

for value_type = 1 : size(values_map_plots, 3)
    subplot(1, 2, value_type);
    
    plot_handle = pcolor(values_map_plots(:, :, value_type));
    %caxis([min(results_comp(:)) max(results_comp(:))]);
    colorbar;
    colormap inferno;
    
    title(measure_labels{plot_measure});
    set(gca, 'Box', 'off');
    set(plot_handle, 'EdgeColor', 'none'); % remove grid outline
end

end
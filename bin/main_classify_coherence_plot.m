%% Load

results_directory = 'workspace_results/';
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat'; % 2.25s x 8 trials coherence
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification_across1.mat'; % 18s x 1 trial coherence

load([results_directory results_filename]);

nNetworks = size(accuracies, 2);
nFlies = size(accuracies, 3);

%% Plot average

figure;
errorbar(...
    frequencies,...
    mean(mean(accuracies, 2), 3),... % Average across channel pairs, then the flies
    0.5*std(mean(accuracies, 2), [], 3) / sqrt(nFlies)); % Standard error across flies
title('Classification (anest/awake) using power');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass 45 100]);

% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);

%% Plot for 1 fly

fly = 11;
figure;
for network = 1 : nNetworks
    plot(frequencies, accuracies(:, network, fly)); hold on;
end
title('Classification (anest/awake) using power');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass 0 100]);

% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);

%% Plot density-normalised values on center/distance mapping

freq_range = (1:42); %(1:83); % corresponding to ~10Hz, check the 'frequencies' vector
freq_range_string = '0-5Hz'; %'0-10Hz';

% Average across frequency range
freq_range_accuracies = permute(mean(accuracies(freq_range, :, :)), [2 3 1]);
freq_range_values = permute(mean(mean(coherencies(freq_range, :, :, :, :), 2), 1), [3 4 5 1 2]); % Average across trials as well

measure_type = 'class'; % 'value' or 'class'
measure_string = 'C';
nChannels = 2;

figure;
set(gcf, 'Position', [0 0 0.3*2100/1.5 300]);
colormap('jet');
channel_sets = networks;

% Get values
if strcmp(measure_type, 'value')
    values = mean(freq_range_values, 3); % Average across flies
    
    % Select 1
    plot_values = values(:, 1) - values(:, 2); fig_tit = 'diff'; cbar_title = [measure_string '_W - ' measure_string '_A']; % Awake - Anest
    %plot_values = values(:, 1)./values(:, 2); fig_tit = 'diff_rel'; cbar_title = [measure_string '_W/' measure_string '_N'];
    %plot_values = values(:, 1); fig_tit = 'wake'; cbar_title = [measure_string '_W']; % Awake
    %plot_values = values(:, 2); fig_tit = 'anest'; cbar_title = [measure_string '_N']; % Anest
    
else % strcmp(measure_type, 'class');
    values = mean(freq_range_accuracies, 2); fig_tit = 'class'; cbar_title = 'class. acc. %';
    %values(values<58) = NaN;
    plot_values = values;
end

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
centers_axis = (min(centers) - (2*centers_delta) : centers_delta : max(centers) + (2*centers_delta)); % - and + delta is for padding
distances_axis = (min(distances) - (2*distances_delta) : distances_delta : max(distances) + (2*distances_delta)); % - and + delta is for padding
[centers_map, distances_map] = meshgrid(centers_axis, distances_axis);

% Create mapped space
values_map = zeros(size(centers_map)); % Will sum all values with the same coordinates
values_map_counter = zeros(size(centers_map)); % Keeps count in each coordinate as to how many values have that coordinate

% Populate mapped space
for value_counter = 1 : length(plot_values)
    x = find(abs(centers_axis - centers(value_counter)) < 0.00001, 1); % This gives the mapped x location
    y = find(distances_axis == distances(value_counter), 1); % This gives the mapped y location
    values_map(y, x) = values_map(y, x) + plot_values(value_counter); % matrix is (rows, columns), corresponding to (y, x)
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
%values_map_plot(values_map_plot < min_value) = NaN; % Remove ghost interpolation (which will affect the colorscale)

%imagesc(values_map_plot); c = colorbar;
plot = pcolor(values_map_plot); c = colorbar;
set(gca, 'color', [0 0 0]); % black background
set(plot, 'EdgeColor', 'none'); % remove grid outline
%axis('xy');

title([num2str(nChannels) 'Ch']);
title(c, cbar_title);
xlabel('set center'); ylabel('set path distance');

figure_name = ['figures_20180606/' fig_tit];
print(figure_name, '-dpng'); % PNG
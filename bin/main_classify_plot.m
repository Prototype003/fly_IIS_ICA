%% Description

%{
For showing classification results
%}

%% General setup

results_directory = 'workspace_results/';
nConditions = 2;

freq_range = (1:42); %(1:83); % corresponding to ~10Hz, check the 'frequencies' vector
freq_range_string = '0-5Hz'; %'0-10Hz';

accuracy_lims = [45 100];

%% Box pot summaries WITHIN
% Mean, median, min, max, for power, coherence, phi-three(2,3,4ch), and phi-star(2,3,4ch), for both within and across flies
% Within - average across channels/sets, variance is across flies

measure_accuracies = []; measure_groups = [];

% Power
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification.mat';
load([results_directory results_filename]);
measure_accuracies = permute(mean(mean(accuracies(freq_range, :, :), 1), 2), [3 1 2]); % average across frequency range and channels
measure_groups = zeros(size(measure_accuracies)) + 1;

% Coherence
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';
load([results_directory results_filename]);
measure_accuracies = [measure_accuracies; permute(mean(mean(accuracies(freq_range, :, :), 1), 2), [3 1 2])]; % average across frequency range and sets
measure_groups = [measure_groups; zeros(size(accuracies, 3), 1) + 2];

% Phi-three
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
load([results_directory results_filename]);
group_counter = 3;
for nChannels = 1 : length(accuracies)
    values = accuracies{nChannels}.accuracies;
    measure_accuracies = [measure_accuracies; permute(mean(values, 1), [2 1])];
    measure_groups = [measure_groups; zeros(size(values, 2), 1) + group_counter];
    group_counter = group_counter + 1;
end

% Phi-star
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
load([results_directory results_filename]);
for nChannels = 1 : length(accuracies)
    values = accuracies{nChannels}.accuracies;
    measure_accuracies = [measure_accuracies; permute(mean(values, 1), [2 1])];
    measure_groups = [measure_groups; zeros(size(values, 2), 1) + group_counter];
    group_counter = group_counter + 1;
end

figure;
boxplot(measure_accuracies, measure_groups);
ylim([49 81]);
ylabel('%');
title('Classification within flies (N=13)');
xticklabels({'1ch P', '2ch C', '2ch Phi3', '3ch Phi3', '4ch Phi3', '2ch Phi*', '3ch Phi*', '4ch Phi*'});

% Plot mean, max, min
hold on;
for measure = 1 : max(measure_groups)
    scatter(measure, mean(measure_accuracies(measure_groups==measure)), 'bx');
    scatter(measure, min(measure_accuracies(measure_groups==measure)), 'b^');
    scatter(measure, max(measure_accuracies(measure_groups==measure)), 'bo');
end

%% Box pot summaries ACROSS
% Mean, median, min, max, for power, coherence, phi-three(2,3,4ch), and phi-star(2,3,4ch), for both within and across flies
% Across - variance is across channels/sets

measure_accuracies = []; measure_groups = [];

% Power
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across1.mat';
load([results_directory results_filename]);
measure_accuracies = permute(mean(accuracies(freq_range, :), 1), [2 1]); % average across frequency range
measure_groups = zeros(size(measure_accuracies)) + 1;

% Coherence
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification_across1.mat';
load([results_directory results_filename]);
measure_accuracies = [measure_accuracies; permute(mean(accuracies(freq_range, :, :), 1), [2 1])]; % average across frequency range and sets
measure_groups = [measure_groups; zeros(size(accuracies, 2), 1) + 2];

% Phi-three
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
load([results_directory results_filename]);
group_counter = 3;
for nChannels = 1 : length(accuracies)
    values = accuracies{nChannels}.accuracies;
    measure_accuracies = [measure_accuracies; values];
    measure_groups = [measure_groups; zeros(size(values, 1), 1) + group_counter];
    group_counter = group_counter + 1;
end

% Phi-star
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
load([results_directory results_filename]);
for nChannels = 1 : length(accuracies)
    values = accuracies{nChannels}.accuracies;
    measure_accuracies = [measure_accuracies; values];
    measure_groups = [measure_groups; zeros(size(values, 1), 1) + group_counter];
    group_counter = group_counter + 1;
end

figure;
boxplot(measure_accuracies, measure_groups);
ylim([20 90]);
ylabel('%');
title('Classification across flies (N=number of channels/sets)');
xticklabels({'1ch P', '2ch C', '2ch Phi3', '3ch Phi3', '4ch Phi3', '2ch Phi*', '3ch Phi*', '4ch Phi*'});

% Plot mean, max, min
hold on;
for measure = 1 : max(measure_groups)
    scatter(measure, mean(measure_accuracies(measure_groups==measure)), 'bx');
    scatter(measure, min(measure_accuracies(measure_groups==measure)), 'b^');
    scatter(measure, max(measure_accuracies(measure_groups==measure)), 'bo');
end

%% Power classification WITHIN FLIES

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification.mat';
load([results_directory results_filename]);

figure;
imagesc(mean(accuracies, 3)); colorbar;
axis('xy');
yticklabels(frequencies([50 100 150 200 250 300 350 400]));
xlabel('channel'); ylabel('frequency bin (range: 0-50Hz)');
title('Mean classification accuracy using power (N=13)');

figure;
errorbar(...
    frequencies,...
    mean(mean(accuracies, 2), 3),... % Average across channels, then the flies
    0.5*std(mean(accuracies, 2), [], 3) / sqrt(size(accuracies, 3))); % standard error across flies
title('Classification (anest/awake) using power (N=13)');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass accuracy_lims]);
% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);

% Get average per channel at set frequency range
figure;
values = permute(mean(mean(accuracies(freq_range, :, :), 1), 3), [2 1 3]);
values_err = permute(std(mean(accuracies(freq_range, :, :), 1), [], 3), [2, 1, 3]) / sqrt(size(accuracies, 3));
%values = permute(mean(mean(mean(powers(freq_range, :, :, 1, :) - powers(freq_range, :, :, 2, :), 1), 2), 5), [3 1 2 4 5]);
errorbar((1:size(accuracies, 2)), values, values_err/2);
xlim([0 16]);
title(['mean accuracy across ' freq_range_string ' (N=13)']);
xlabel('channel'); ylabel('%');

%% Power classification ACROSS FLIES

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across1.mat';
load([results_directory results_filename]);

figure;
imagesc(mean(accuracies, 3)); colorbar;
axis('xy');
yticklabels(frequencies([50 100 150 200 250 300 350 400]));
xlabel('channel'); ylabel('frequency bin (range: 0-50Hz)');
title('Mean classification accuracy using power');

figure;
errorbar(...
    frequencies,...
    mean(accuracies, 2),... % Average across channels
    0.5*std(mean(accuracies, 2), [], 3) / sqrt(size(accuracies, 2))); % standard error across flies
title('Classification (anest/awake) using power (N=15 channels)');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass accuracy_lims]);
% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);

% Get average per channel at set frequency range
figure;
values = permute(mean(mean(accuracies(freq_range, :, :), 1), 3), [2 1 3]);
values_err = permute(std(mean(accuracies(freq_range, :, :), 1), [], 3), [2, 1, 3]) / sqrt(size(accuracies, 3));
errorbar((1:size(accuracies, 2)), values, values_err/2);
xlim([0 16]);
title(['mean accuracy across ' freq_range_string]);
xlabel('channel'); ylabel('%');

%% Plot density-normalised values on center/distance mapping

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';
load([results_directory results_filename]);

figure;
set(gcf, 'Position', [0 0 0.3*2100/1.5 300]);
colormap('jet');
for nChannels_counter = 1 : length(phis)
    nChannels = phis{nChannels_counter}.nChannels;
    
    channel_sets = phis{nChannels_counter}.channel_sets;
    
    % Get values
    if strcmp(measure_type, 'value')
        values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, flies, :, tau), 2), 3), [1 4 2 3]);
        
        % Select 1
        plot_values = values(:, 1) - values(:, 2); fig_tit = 'diff'; cbar_title = [measure_string '_W - ' measure_string '_A']; % Awake - Anest
        %plot_values = values(:, 1)./values(:, 2); fig_tit = 'diff_rel'; cbar_title = [measure_string '_W/' measure_string '_N'];
        %plot_values = values(:, 1); fig_tit = 'wake'; cbar_title = [measure_string '_W']; % Awake
        %plot_values = values(:, 2); fig_tit = 'anest'; cbar_title = [measure_string '_N']; % Anest
        
    else % strcmp(measure_type, 'class');
        values = mean(accuracies{nChannels_counter}.accuracies, 2); fig_tit = 'class'; cbar_title = 'class. acc. %';
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
    
    subplot(1, length(phis), nChannels_counter);
    %imagesc(values_map_plot); c = colorbar;
    plot = pcolor(values_map_plot); c = colorbar;
    set(gca, 'color', [0 0 0]); % black background
    set(plot, 'EdgeColor', 'none'); % remove grid outline
    %axis('xy');
    
    title([num2str(nChannels) 'Ch']);
    title(c, cbar_title);
    xlabel('set center'); ylabel('set path distance');
end

figure_name = ['figures_20180606/' fig_tit];
print(figure_name, '-dpng'); % PNG

%% Coherence classification WITHIN FLIES

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';
load([results_directory results_filename]);

figure;
imagesc(mean(accuracies, 3)); colorbar;
axis('xy');
yticklabels(frequencies([50 100 150 200 250 300 350 400]));
xlabel('channel pair'); ylabel('frequency bin (range: 0-50Hz)');
title('Mean classification accuracy using coherency (N=13)');

figure;
errorbar(...
    frequencies,...
    mean(mean(accuracies, 2), 3),... % Average across channel pairs, then the flies
    0.5*std(mean(accuracies, 2), [], 3) / sqrt(size(accuracies, 3))... % standard error across flies
    );
title('Classification (anest/awake) using coherence (N=13)');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass accuracy_lims(1) accuracy_lims(2)]);
line(chronux_params.fpass, [100/nConditions 100/nConditions]); % Chance level

% % Accuracy as a function of global channel distance
% figure;
% set_distances = channel_set_distances(networks);
% for f = 1 : size(accuracies)
%     scatter(set_distances, mean(accuracies(f, :, :), 3), '.'); hold on;
% end

% % Plot actual values to see rough separation between conditions
% figure;
% for fly = 1 : 1%size(coherencies, 5)
%     for trial = 1 : size(coherencies, 2)
%         scatter((1:410), coherencies(:, trial, 70, 1, fly), 'r.'); hold on;
%         scatter((1:410), coherencies(:, trial, 70, 2, fly), 'b.');
%     end
% end

%% Coherence classification WITHIN FLIES bar graph collapse

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';
load([results_directory results_filename]);

% % Collapse at all frequencies
% channel_sets = networks;
% values_collapsed = zeros(max(channel_sets(:)), size(accuracies, 3), size(accuracies, 1));
% for fly = 1 : size(accuracies, 3)
%     
%     for frequency = 1 : size(accuracies, 1)
%         
%         values = permute(accuracies(frequency, :, fly), [2 1 3]);
%         
%         % Sum accuracies for each channel (sum across networks which contain the channel)
%         set_counters = zeros(max(channel_sets(:)), 1);
%         for channel = 1 : max(channel_sets(:))
%             for channel_set = 1 : size(channel_sets, 1)
%                 if any(channel_sets(channel_set, :) == channel)
%                     set_counters(channel) = set_counters(channel) + 1;
%                     values_collapsed(channel, fly, frequency) = values_collapsed(channel, fly, frequency) + values(channel_set);
%                 end
%             end
%         end
%         
%         % Average accuracies
%         values_collapsed(:, fly, frequency) = values_collapsed(:, fly, frequency) ./ set_counters;
%         
%     end
%     
% end

% Collapse after averaging across frequency range
channel_sets = networks;
freq_range_values = permute(mean(accuracies(freq_range, :, :)), [2 3 1]);
channel_values = zeros(max(channel_sets(:)), size(accuracies, 3));
set_counters = zeros(size(channel_values));
for channel = 1 : size(channel_values, 1)
    for channel_set = 1 : size(channel_sets, 1)
        if any(channel_sets(channel_set, :) == channel)
            set_counters(channel, :) = set_counters(channel) + 1;
            channel_values(channel, :) = channel_values(channel, :) + freq_range_values(channel_set, :);
        end
    end
end
% Average accuracies
channel_values = channel_values ./ set_counters;

% Plot average per channel (averaged across flies)
figure;
values = mean(channel_values, 2);
values_err = std(channel_values, [], 2) / sqrt(size(channel_values, 2));
errorbar((1:size(channel_values, 1)), values, values_err/2);
xlim([0 16]);
title(['mean accuracy across ' freq_range_string]);
xlabel('channel'); ylabel('%');

% Path distance plot
figure; colormap('jet');
centers = mean(channel_sets, 2); % Mean across channels in each set
distances = channel_set_distances(channel_sets);
centers_offset_rand = (-0.1+(0.1--0.1).*rand(length(centers), 1)); % Values can vary by just under +- half of the minimum difference in centers
distances_offset_rand = (-0.4+(0.4--0.4).*rand(length(distances), 1)); % Values can vary by just under +- half of the minimum difference in distances
scatter(centers+centers_offset_rand, distances+distances_offset_rand, 100, mean(freq_range_values, 2), nChannels_shapes(nChannels_counter)); c = colorbar;hold on;
ax = gca; ax.Color = 'k';
xlabel('set center');
ylabel('set path distance');
title(c, 'class. acc. %');

%% Coherence classification ACROSS FLIES

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification_across1.mat';
load([results_directory results_filename]);

figure;
imagesc(mean(accuracies, 3)); colorbar;
axis('xy');
yticklabels(frequencies([50 100 150 200 250 300 350 400]));
xlabel('channel pair'); ylabel('frequency bin (range: 0-50Hz)');
title('Mean classification accuracy using coherency (N=13)');

figure;
errorbar(...
    frequencies,...
    mean(accuracies, 2),... % Average across channel pairs
    0.5*std(accuracies, [], 2) / sqrt(size(accuracies, 2))... % standard error across channel sets
    );
title('Classification (anest/awake) using coherence (N=105 channel pairs)');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass accuracy_lims(1) accuracy_lims(2)]);
line(chronux_params.fpass, [100/nConditions 100/nConditions]); % Chance level

% Collapse after averaging across frequency range
channel_sets = networks;
freq_range_values = permute(mean(accuracies(freq_range, :, :)), [2 3 1]);
channel_values = zeros(max(channel_sets(:)), size(accuracies, 3));
set_counters = zeros(size(channel_values));
for channel = 1 : size(channel_values, 1)
    for channel_set = 1 : size(channel_sets, 1)
        if any(channel_sets(channel_set, :) == channel)
            set_counters(channel, :) = set_counters(channel) + 1;
            channel_values(channel, :) = channel_values(channel, :) + freq_range_values(channel_set, :);
        end
    end
end
% Average accuracies
channel_values = channel_values ./ set_counters;

% Plot average per channel (averaged across flies)
figure;
values = mean(channel_values, 2);
values_err = std(channel_values, [], 2) / sqrt(size(channel_values, 2));
errorbar((1:size(channel_values, 1)), values, values_err/2);
xlim([0 16]);
title(['mean accuracy across ' freq_range_string]);
xlabel('channel'); ylabel('%');

% Path distance plot
figure; colormap('jet');
centers = mean(channel_sets, 2); % Mean across channels in each set
distances = channel_set_distances(channel_sets);
centers_offset_rand = (-0.1+(0.1--0.1).*rand(length(centers), 1)); % Values can vary by just under +- half of the minimum difference in centers
distances_offset_rand = (-0.4+(0.4--0.4).*rand(length(distances), 1)); % Values can vary by just under +- half of the minimum difference in distances
scatter(centers+centers_offset_rand, distances+distances_offset_rand, 100, mean(freq_range_values, 2), nChannels_shapes(nChannels_counter)); c = colorbar;hold on;
ax = gca; ax.Color = 'k';
xlabel('set center');
ylabel('set path distance');
title(c, 'class. acc. %');

%% Coherence classification ACROSS FLIES bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average classification accuracy across channel sets which include the channel

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification_across1.mat';
load([results_directory results_filename]);

colours = 'rgb';
channel_sets = networks;
values = accuracies;

values_collapsed = zeros(size(accuracies, 1), max(channel_sets(:)));
for frequency = 1 : size(accuracies, 1)
    
    % Sum accuracies for each channel (sum across networks which contain the channel)
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed(frequency, channel) = values_collapsed(frequency, channel) + values(frequency, channel_set);
            end
        end
    end
    
    % Average accuracies
    values_collapsed(frequency, :) = values_collapsed(frequency, :) ./ set_counters;
    
end

% Plot
figure;
imagesc(values_collapsed); colorbar;
axis('xy');

%% Phi-three classification WITHIN FLIES

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
load([results_directory results_filename]);

figure;
for nChannels_counter = 1 : length(accuracies)
    subplot(1, length(accuracies), nChannels_counter);
    imagesc(accuracies{nChannels_counter}.accuracies'); colorbar;
    title(['Classification (anest/awake) ' num2str(nChannels_counter+1) 'ch']);
    xlabel('channel set');
    ylabel('fly');
end

figure;
for nChannels_counter = 1 : length(accuracies)
    subplot(1, length(accuracies), nChannels_counter);
%     errorbar(...
%         (1:size(accuracies{nChannels_counter}.accuracies, 1)),...
%         mean(accuracies{nChannels_counter}.accuracies, 2),...
%         0.5*std(accuracies{nChannels_counter}.accuracies, [], 2) / sqrt(size(accuracies{nChannels_counter}.accuracies, 2))...
%         );
    plot(...
        (1:size(accuracies{nChannels_counter}.accuracies, 1)),...
        mean(accuracies{nChannels_counter}.accuracies, 2)...
        );
    title(['Phi3 Classification ' num2str(nChannels_counter+1) 'ch']);
    xlabel('channel set');
    ylabel('% correct');
    axis([1 size(accuracies{nChannels_counter}.accuracies, 1) accuracy_lims]);
    % Chance level
    line([1 size(accuracies{nChannels_counter}.accuracies, 1)], [100/nConditions 100/nConditions]);
end

% Accuracy as a function of global channel distance (per fly)
nChannels_counter = 3;
set_distances = channel_set_distances(accuracies{nChannels_counter}.channel_sets);
figure;
for fly = 1:13
    subplot(4, 4, fly);
    plot(set_distances, accuracies{nChannels_counter}.accuracies(:, fly), '.'); hold on;
    title(['fly ' num2str(fly)]);
    axis([0 55 0 100]);
end

% Accuracy as a function of global channel distance
figure;
for nChannels_counter = 1 : length(accuracies)
    subplot(1, length(accuracies), nChannels_counter);
    set_distances = channel_set_distances(accuracies{nChannels_counter}.channel_sets);
    scatter(set_distances, mean(accuracies{nChannels_counter}.accuracies(:, :), 2), '.'); hold on; % Mean accuracy across flies
    [a, b] = corr(set_distances, mean(accuracies{nChannels_counter}.accuracies(:, :), 2))
    axis([0 55 50 100]);
    xlabel('sum of channel distances');
    ylabel('accuracy %');
end

%% Phi-three classification WITHIN FLIES 2D (2-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
load([results_directory results_filename]);

% Build visual matrix
channel_sets = accuracies{1}.channel_sets + 1; % +1 is to convert from 0-indexing to 1-indexing
value_map = zeros(max(channel_sets(:)));
for value = 1 : length(accuracies{1}.accuracies)
    value_map(channel_sets(value, 1), channel_sets(value, 2)) = mean(accuracies{1}.accuracies(value, :));
end
% Plot visual matrix
figure;
imagesc(value_map, [min(value_map(:)) max(value_map(:))]); cbar = colorbar; ylabel(cbar, 'accuracy');
xlabel('channel'); ylabel('channel');

% Or plot value by value
% Colour scaling from accuracy to colour: https://math.stackexchange.com/questions/914823/shift-numbers-into-a-different-range
% But matlab linearly scales vectors to the current colormap anyway
figure;
colormap('jet');
values = mean(accuracies{1}.accuracies, 2);
values(values<58) = NaN;
scatter(channel_sets(:, 1), channel_sets(:, 2), values, values, 'filled'); hold on;
xlabel('channel'); ylabel('channel'); title('classification accuracy 2ch');
colorbar;
%caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
caxis([min(values) max(values)]);
axis([0 16 0 16]);

%% Phi-three classifiation WITHIN FLIES 3D (3-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
load([results_directory results_filename]);

channel_sets = accuracies{2}.channel_sets+1;

figure;
colormap('jet');
values = mean(accuracies{2}.accuracies, 2);
values(values<58) = NaN;
scatter3(channel_sets(:, 1), channel_sets(:, 2), channel_sets(:, 3), values/10, values, 'filled');
xlabel('channel'); ylabel('channel'); zlabel('channel'); title('classification accuracy 3ch');
axis([0 16 0 16 0 16]);
%caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
colorbar;
caxis([min(values) max(values)]);
view([4 25]); % [azimuth elevation]

%% Phi-three classification WITHIN FLIES 4D (4-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
load([results_directory results_filename]);

values = mean(accuracies{3}.accuracies, 2);
channel_sets = double(accuracies{3}.channel_sets);

figure;
colormap('jet');
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices);
    values(values<58) = NaN;
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), values_sliced, values_sliced, '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    if fourth_d == min(channels_4d)
        colorbar;
    end
    axis([0 16 0 16 0 16]);
    %caxis([58 max(accuracies{3}.accuracies)]);
    %caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
    caxis([min(values) max(values)]);
    view([6 25]); % [azimuth elevation]
end

%% Phi-three classification WITHIN FLIES bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average classification accuracy across channel sets which include the channel

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
load([results_directory results_filename]);

python_indexing = [1 1 0];

figure;
colours = 'rgb';
for nChannels_counter = 1 : length(accuracies)
    nChannels = accuracies{nChannels_counter}.nChannels;
    channel_sets = double(accuracies{nChannels_counter}.channel_sets) + python_indexing(nChannels_counter);
    values = mean(accuracies{nChannels_counter}.accuracies, 2);
    
    % Sum accuracies for each channel (sum across networks which contain the channel)
    values_collapsed = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
            end
        end
    end
    
    % Average accuracies
    values_collapsed = values_collapsed ./ set_counters;
    % Standard error
    values_collapsed_deviations = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
            end
        end
    end
    values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
    values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;
    
    % Plot
    %subplot(1, length(accuracies), nChannels_counter);
    plot((1:length(values_collapsed)), values_collapsed, colours(nChannels_counter)); hold on;
    %bar((1:length(values_collapsed)), values_collapsed);
    %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr);
    axis([0 16 50 75]);
    legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
    ylabel('%'); xlabel('channel'); title('Mean accuracy across sets containing channel X');
end

%% Phi-three classification ACROSS FLIES

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
load([results_directory results_filename]);

figure;
for nChannels_counter = 1 : length(accuracies)
    subplot(1, length(accuracies), nChannels_counter);
    plot((1:length(accuracies{nChannels_counter}.accuracies)), accuracies{nChannels_counter}.accuracies); hold on;
    line([1 length(accuracies{nChannels_counter}.accuracies)], [100/nConditions 100/nConditions]);
    axis([1 length(accuracies{nChannels_counter}.accuracies) 20 100]);
    title(['Phi3 Classification ' num2str(nChannels_counter+1) 'ch']);
    xlabel('channel set');
    ylabel('% correct');
end

%% Phi-three classification ACROSS FLIES 2D (2-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
load([results_directory results_filename]);

% Build visual matrix
channel_sets = accuracies{1}.channel_sets + 1; % +1 is to convert from 0-indexing to 1-indexing
value_map = zeros(max(channel_sets(:)));
for value = 1 : length(accuracies{1}.accuracies)
    value_map(channel_sets(value, 1), channel_sets(value, 2)) = accuracies{1}.accuracies(value);
end
% Plot visual matrix
figure;
imagesc(value_map, [min(accuracies{1}.accuracies) max(accuracies{1}.accuracies)]); cbar = colorbar; ylabel(cbar, 'accuracy');
xlabel('channel'); ylabel('channel');

% Or plot value by value
% Colour scaling from accuracy to colour: https://math.stackexchange.com/questions/914823/shift-numbers-into-a-different-range
% But matlab linearly scales vectors to the current colormap anyway
figure;
colormap('jet');
values = (accuracies{1}.accuracies);
values(values<58) = NaN;
scatter(channel_sets(:, 1), channel_sets(:, 2), values, values, 'filled'); hold on;
xlabel('channel'); ylabel('channel'); title('classification accuracy 2ch');
colorbar;
%caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
caxis([min(values) max(values)]);
axis([0 16 0 16]);

%% Phi-three classifiation ACROSS FLIES 3D (3-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
load([results_directory results_filename]);

channel_sets = accuracies{2}.channel_sets+1;

figure;
colormap('jet');
values = (accuracies{2}.accuracies);
values(values<58) = NaN;
scatter3(channel_sets(:, 1), channel_sets(:, 2), channel_sets(:, 3), values/10, values, 'filled');
xlabel('channel'); ylabel('channel'); zlabel('channel'); title('classification accuracy 3ch');
axis([0 16 0 16 0 16]);
caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
colorbar;
caxis([min(values) max(values)]);
view([4 25]); % [azimuth elevation]

%% Phi-three classification ACROSS FLIES 4D (4-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
load([results_directory results_filename]);

values = (accuracies{3}.accuracies);
channel_sets = double(accuracies{3}.channel_sets)+1;

figure;
colormap('jet');
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices);
    values(values<58) = NaN;
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), values_sliced, values_sliced, '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    if fourth_d == min(channels_4d)
        colorbar;
    end
    axis([0 16 0 16 0 16]);
    caxis([58 max(accuracies{3}.accuracies)]);
    %caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
    %caxis([min(values) max(values)]);
    view([6 25]); % [azimuth elevation]
end

%% Phi-three classification ACROSS FLIES bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average classification accuracy across channel sets which include the channel

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
load([results_directory results_filename]);

figure;
colours = 'rgb';
for nChannels_counter = 1 : length(accuracies)
    nChannels = accuracies{nChannels_counter}.nChannels;
    channel_sets = double(accuracies{nChannels_counter}.channel_sets) + 1;
    values = accuracies{nChannels_counter}.accuracies;
    
    % Sum accuracies for each channel (sum across networks which contain the channel)
    values_collapsed = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
            end
        end
    end
    
    % Average accuracies
    values_collapsed = values_collapsed ./ set_counters;
    % Standard error
    values_collapsed_deviations = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
            end
        end
    end
    values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
    values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;
    
    % Plot
    %subplot(1, length(accuracies), nChannels_counter);
    plot((1:length(values_collapsed)), values_collapsed, colours(nChannels_counter)); hold on;
    %bar((1:length(values_collapsed)), values_collapsed);
    %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr);
    axis([0 16 50 75]);
    legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
    ylabel('%'); xlabel('channel'); title('Mean accuracy across sets containing channel X');
end

%% LOAD Phi-three values GLOBAL/NON-GLOBAL FLIES 2D, 3D, and 4D plots LOAD

% Global TPM (18s) was used in across fly analysis
% Non-global TPM (2.25s per trial) was used in within fly analysis
% All classification analyses were conducted at tau = 4ms

global_tpm = 0; % 0: 2.25s TPMs; 1: 18s TPMs

if global_tpm == 1
    % For across fly analysis we want one value per fly, so we use 18s TPM
    disp('loading');
    data_directory = 'results/';
    data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree.mat'];
    load([data_directory data_filename]);
    disp('loaded');
else % global_tpm == 0
    % For within fly analysis we use 16 values per fly (8 awake and 8 anest), so we used 2.25s TPMs
    phis = cell(1, 3);
    data_directory = 'results/';
    data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t2_phithree_nonGlobal.mat'];
    disp('loading 2ch');
    tmp = load([data_directory data_filename]);
    phis{1} = tmp.phis{1};
    data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal_tau4.mat'];
    disp('loading 3ch');
    tmp = load([data_directory data_filename]);
    phis{2} = tmp.phis{1};
    data_directory = 'results/';
    data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal_tau4.mat'];
    disp('loading 4ch');
    tmp = load([data_directory data_filename]);
    phis{3} = tmp;
    phis{3}.channel_sets = phis{3}.phis{1}.channel_sets-1; % Make indexing consistent with other nChannels (python indexing)
    clear tmp
    disp('loaded');
end
tau = 1;

%% Phi-three values ACROSS FLIES 2D plot

% Build visual matrix
channel_sets = phis{1}.channel_sets + 1;
value_map = zeros(max(channel_sets(:)));
values = permute(mean(mean(phis{1}.phi_threes(:, :, :, :, tau), 2), 3), [1 4 2 3]);
for value = 1 : size(values, 1)
    value_map(channel_sets(value, 1), channel_sets(value, 2)) = values(value, 1); % - values(value, 2)
end
figure;
colormap('jet');
imagesc(value_map); cbar = colorbar; ylabel(cbar, '\Phi');
xlabel('channel'); ylabel('channel');
title('Wake');

%% Phi-three values ACROSS FLIES 3D plot

values = permute(mean(mean(phis{2}.phi_threes(:, :, :, :, tau), 2), 3), [1 4 2 3]);
channel_sets = phis{2}.channel_sets + 1;
figure; colormap('jet');
scatter3(channel_sets(:, 1), channel_sets(:, 2), channel_sets(:, 3), 75, values(:, 1), '.');
cbar = colorbar; ylabel(cbar, '\Phi');
xlabel('channel'); ylabel('channel'); zlabel('channel');
axis([0 16 0 16 0 16]);
view([4 25]); % [azimuth elevation]

%% Phi-three values ACROSS FLIES 4D plot

values = permute(mean(mean(phis{3}.phi_threes(:, :, :, :, tau), 2), 3), [1 4 2 3]);
channel_sets = double(phis{3}.channel_sets) + 1;

% AWAKE - ANEST
figure;
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices, :);
    
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), 75, values_sliced(:, 1) - values_sliced(:, 2), '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    axis([0 16 0 16 0 16]);
    caxis([min(min(values(:, 1) - values(:, 2))) max(max(values(:, 1) - values(:, 2)))]);
    view([6 25]); % [azimuth elevation]
end

% AWAKE
figure; colormap('jet');
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices, :);
    
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), 75, values_sliced(:, 1), '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    axis([0 16 0 16 0 16]);
    %caxis([min(values(:, 2)) max(values(:, 1))]);
    view([6 25]); % [azimuth elevation]
end

% ANEST
figure;
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices, :);
    
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), 75, values_sliced(:, 2), '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    axis([0 16 0 16 0 16]);
    caxis([min(values(:, 2)) max(values(:, 1))]);
    view([6 25]); % [azimuth elevation]
end

%% Phi-three values ACROSS FLIES bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average phi value across channel sets which include the channel

figure;
colours = 'rgb';
for nChannels_counter = 1 : length(phis)
    nChannels = phis{nChannels_counter}.nChannels;
    channel_sets = double(phis{nChannels_counter}.channel_sets) + 1;
    values = permute(mean(mean(phis{nChannels_counter}.phi_threes(:, :, :, :, tau), 2), 3), [1 4 2 3]);
    values = values(:, 2);% - values(:, 2);
    
    % Sum accuracies for each channel (sum across networks which contain the channel)
    values_collapsed = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
            end
        end
    end
    
    % Average phis
    values_collapsed = values_collapsed ./ set_counters;
    % Standard error
    values_collapsed_deviations = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
            end
        end
    end
    values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
    values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;
    
    % Plot
    %subplot(1, length(accuracies), nChannels_counter);
    plot((1:length(values_collapsed)), values_collapsed, colours(nChannels_counter)); hold on;
    %bar((1:length(values_collapsed)), values_collapsed);
    %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr, colours(nChannels_counter)); hold on;
    axis([0 16 0 0.03]);
    legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
    ylabel('\Phi'); xlabel('channel'); title('Mean phi3 across sets containing channel X');
end

%% Phi-three values WITHIN FLIES bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average phi value across channel sets which include the channel

figure;
colours = 'rgb';
for fly = 1 : 13
    subplot(4, 4, fly);
    for nChannels_counter = 1 : length(phis)
        nChannels = phis{nChannels_counter}.nChannels;
        channel_sets = double(phis{nChannels_counter}.channel_sets) + 1;
        values = permute(mean(mean(phis{nChannels_counter}.phi_threes(:, :, fly, :, tau), 2), 3), [1 4 2 3]);
        values = values(:, 2);% - values(:, 2);
        
        % Sum accuracies for each channel (sum across networks which contain the channel)
        values_collapsed = zeros(1, max(channel_sets(:)));
        set_counters = zeros(1, max(channel_sets(:)));
        for channel = 1 : max(channel_sets(:))
            for channel_set = 1 : size(channel_sets, 1)
                if any(channel_sets(channel_set, :) == channel)
                    set_counters(channel) = set_counters(channel) + 1;
                    values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
                end
            end
        end
        
        % Average phis
        values_collapsed = values_collapsed ./ set_counters;
        % Standard error
        values_collapsed_deviations = zeros(1, max(channel_sets(:)));
        set_counters = zeros(1, max(channel_sets(:)));
        for channel = 1 : max(channel_sets(:))
            for channel_set = 1 : size(channel_sets, 1)
                if any(channel_sets(channel_set, :) == channel)
                    set_counters(channel) = set_counters(channel) + 1;
                    values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
                end
            end
        end
        values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
        values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;
        
        % Plot
        %subplot(1, length(accuracies), nChannels_counter);
        plot((1:length(values_collapsed)), values_collapsed, colours(nChannels_counter)); hold on;
        %bar((1:length(values_collapsed)), values_collapsed);
        %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr, colours(nChannels_counter)); hold on;
        axis([0 16 0 0.05]);
        legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
        ylabel('\Phi'); xlabel('channel'); title('Mean phi3 across sets containing channel X');
    end
end

%% Phi-three values ACROSS FLIES 4ch vs mean of channel means

% Get average 2-channel phi for each channel
nChannels_counter = 2;
nChannels = phis{nChannels_counter}.nChannels;
channel_sets = double(phis{nChannels_counter}.channel_sets) + 1;
values = permute(mean(mean(phis{nChannels_counter}.phi_threes(:, :, :, :, tau), 2), 3), [1 4 2 3]);
values = values(:, 1) - values(:, 2); % WHY?????

% Sum accuracies for each channel (sum across networks which contain the channel)
values_collapsed = zeros(1, max(channel_sets(:)));
set_counters = zeros(1, max(channel_sets(:)));
for channel = 1 : max(channel_sets(:))
    for channel_set = 1 : size(channel_sets, 1)
        if any(channel_sets(channel_set, :) == channel)
            set_counters(channel) = set_counters(channel) + 1;
            values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
        end
    end
end

% Average phis
values_collapsed = values_collapsed ./ set_counters;
% Standard error
values_collapsed_deviations = zeros(1, max(channel_sets(:)));
set_counters = zeros(1, max(channel_sets(:)));
for channel = 1 : max(channel_sets(:))
    for channel_set = 1 : size(channel_sets, 1)
        if any(channel_sets(channel_set, :) == channel)
            set_counters(channel) = set_counters(channel) + 1;
            values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
        end
    end
end
values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;

figure;
colours = 'rb';
subplot_counter = 1;
for nChannels_counter = 2 : length(phis)
    subplot(1, 2, nChannels_counter - 1);
    for condition = 1 : 2
        nChannels = phis{nChannels_counter}.nChannels;
        channel_sets = double(phis{nChannels_counter}.channel_sets) + 1;
        values = mean(mean(phis{nChannels_counter}.phi_threes(:, :, :, condition, tau), 2), 3);
        
        % Find mean of mean 2-channel phis for each network
        values_2ch_mean = zeros(size(values));
        for channel_set = 1 : size(channel_sets, 1)
            channels = channel_sets(channel_set, :);
            values_2ch_mean(channel_set) = mean(values_collapsed(channels));
        end
        
        % Plot
        scatter(values_2ch_mean, values, 100, colours(condition), '.'); hold on;
        disp([num2str(nChannels_counter+1) 'Ch, condition ' num2str(condition)]);
        [r, p] = corr(values_2ch_mean, values)
    end
    title([num2str(nChannels) ' channels']);
    xlabel('mean 2-ch phi3 across channels in set');
    ylabel('phi3 value');
    legend('awake', 'anest', 'Location', 'northeastoutside');
end

%% Phi-star classification WITHIN FLIES

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
load([results_directory results_filename]);

figure;
for nChannels_counter = 1 : length(accuracies)
    subplot(1, length(accuracies), nChannels_counter);
    imagesc(accuracies{nChannels_counter}.accuracies'); colorbar;
    title(['Classification (anest/awake) ' num2str(nChannels_counter+1) 'ch']);
    xlabel('channel set');
    ylabel('fly');
end

figure;
for nChannels_counter = 1 : length(accuracies)
    subplot(1, length(accuracies), nChannels_counter);
    %     errorbar(...
    %         (1:size(accuracies{nChannels_counter}.accuracies, 1)),...
    %         mean(accuracies{nChannels_counter}.accuracies, 2),...
    %         0.5*std(accuracies{nChannels_counter}.accuracies, [], 2) / sqrt(size(accuracies{nChannels_counter}.accuracies, 2))...
    %         );
    plot(...
        (1:size(accuracies{nChannels_counter}.accuracies, 1)),...
        mean(accuracies{nChannels_counter}.accuracies, 2)...
        );
    title(['Phi* Classification ' num2str(nChannels_counter+1) 'ch']);
    xlabel('channel set');
    ylabel('% correct');
    axis([1 size(accuracies{nChannels_counter}.accuracies, 1) accuracy_lims]);
    % Chance level
    line([1 size(accuracies{nChannels_counter}.accuracies, 1)], [100/nConditions 100/nConditions]);
end

%% Phi-star classification WITHIN FLIES 2D (2-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
load([results_directory results_filename]);

% Build visual matrix
channel_sets = accuracies{1}.channel_sets;
value_map = zeros(max(channel_sets(:)));
for value = 1 : length(accuracies{1}.accuracies)
    value_map(channel_sets(value, 1), channel_sets(value, 2)) = mean(accuracies{1}.accuracies(value, :));
end
% Plot visual matrix
figure;
imagesc(value_map, [min(value_map(:)) max(value_map(:))]); cbar = colorbar; ylabel(cbar, 'accuracy');
xlabel('channel'); ylabel('channel');

% Or plot value by value
% Colour scaling from accuracy to colour: https://math.stackexchange.com/questions/914823/shift-numbers-into-a-different-range
% But matlab linearly scales vectors to the current colormap anyway
figure;
colormap('jet');
values = mean(accuracies{1}.accuracies, 2);
values(values<58) = NaN;
scatter(channel_sets(:, 1), channel_sets(:, 2), values, values, 'filled'); hold on;
xlabel('channel'); ylabel('channel'); title('classification accuracy 2ch');
colorbar;
%caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
caxis([min(values) max(values)]);
axis([0 16 0 16]);

%% Phi-star classifiation WITHIN FLIES 3D (3-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
load([results_directory results_filename]);

channel_sets = accuracies{2}.channel_sets;

figure;
colormap('jet');
values = mean(accuracies{2}.accuracies, 2);
values(values<58) = NaN;
scatter3(channel_sets(:, 1), channel_sets(:, 2), channel_sets(:, 3), values/10, values, 'filled');
xlabel('channel'); ylabel('channel'); zlabel('channel'); title('classification accuracy 3ch');
axis([0 16 0 16 0 16]);
%caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
colorbar;
caxis([min(values) max(values)]);
view([4 25]); % [azimuth elevation]

%% Phi-star classification WITHIN FLIES 4D (4-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
load([results_directory results_filename]);

values = mean(accuracies{3}.accuracies, 2);
channel_sets = double(accuracies{3}.channel_sets);

figure;
colormap('jet');
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices);
    values(values<58) = NaN;
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), values_sliced, values_sliced, '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    if fourth_d == min(channels_4d)
        colorbar;
    end
    axis([0 16 0 16 0 16]);
    %caxis([58 max(accuracies{3}.accuracies)]);
    %caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
    caxis([min(values) max(values)]);
    view([6 25]); % [azimuth elevation]
end

%% Phi-star classification WITHIN FLIES bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average classification accuracy across channel sets which include the channel

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
load([results_directory results_filename]);

figure;
colours = 'rgb';
for nChannels_counter = 1 : length(accuracies)
    nChannels = accuracies{nChannels_counter}.nChannels;
    channel_sets = double(accuracies{nChannels_counter}.channel_sets);
    values = mean(accuracies{nChannels_counter}.accuracies, 2);
    
    % Sum accuracies for each channel (sum across networks which contain the channel)
    values_collapsed = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
            end
        end
    end
    
    % Average accuracies
    values_collapsed = values_collapsed ./ set_counters;
    % Standard error
    values_collapsed_deviations = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
            end
        end
    end
    values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
    values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;
    
    % Plot
    %subplot(1, length(accuracies), nChannels_counter);
    plot((1:length(values_collapsed)), values_collapsed, colours(nChannels_counter)); hold on;
    %bar((1:length(values_collapsed)), values_collapsed);
    %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr);
    axis([0 16 50 75]);
    legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
    ylabel('%'); xlabel('channel'); title('Mean accuracy across sets containing channel X');
end

%% Phi-star classification ACROSS FLIES

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
load([results_directory results_filename]);

figure;
for nChannels_counter = 1 : length(accuracies)
    subplot(1, length(accuracies), nChannels_counter);
    imagesc(accuracies{nChannels_counter}.accuracies'); colorbar;
    title(['Classification (anest/awake) ' num2str(nChannels_counter+1) 'ch']);
    xlabel('channel set');
    ylabel('fly');
end

figure;
for nChannels_counter = 1 : length(accuracies)
    subplot(1, length(accuracies), nChannels_counter);
    %     errorbar(...
    %         (1:size(accuracies{nChannels_counter}.accuracies, 1)),...
    %         mean(accuracies{nChannels_counter}.accuracies, 2),...
    %         0.5*std(accuracies{nChannels_counter}.accuracies, [], 2) / sqrt(size(accuracies{nChannels_counter}.accuracies, 2))...
    %         );
    plot(...
        (1:size(accuracies{nChannels_counter}.accuracies, 1)),...
        mean(accuracies{nChannels_counter}.accuracies, 2)...
        );
    title(['Phi* Classification ' num2str(nChannels_counter+1) 'ch']);
    xlabel('channel set');
    ylabel('% correct');
    axis([1 size(accuracies{nChannels_counter}.accuracies, 1) accuracy_lims]);
    % Chance level
    line([1 size(accuracies{nChannels_counter}.accuracies, 1)], [100/nConditions 100/nConditions]);
end

%% Phi-star classification ACROSS FLIES 2D (2-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
load([results_directory results_filename]);

% Build visual matrix
channel_sets = accuracies{1}.channel_sets;
value_map = zeros(max(channel_sets(:)));
for value = 1 : length(accuracies{1}.accuracies)
    value_map(channel_sets(value, 1), channel_sets(value, 2)) = accuracies{1}.accuracies(value);
end
% Plot visual matrix
figure;
imagesc(value_map, [min(accuracies{1}.accuracies) max(accuracies{1}.accuracies)]); cbar = colorbar; ylabel(cbar, 'accuracy');
xlabel('channel'); ylabel('channel');

% Or plot value by value
% Colour scaling from accuracy to colour: https://math.stackexchange.com/questions/914823/shift-numbers-into-a-different-range
% But matlab linearly scales vectors to the current colormap anyway
figure;
colormap('jet');
values = (accuracies{1}.accuracies);
values(values<58) = NaN;
scatter(channel_sets(:, 1), channel_sets(:, 2), values, values, 'filled'); hold on;
xlabel('channel'); ylabel('channel'); title('classification accuracy 2ch');
colorbar;
%caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
caxis([min(values) max(values)]);
axis([0 16 0 16]);

%% Phi-star classifiation ACROSS FLIES 3D (3-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
load([results_directory results_filename]);

channel_sets = accuracies{2}.channel_sets;

figure;
colormap('jet');
values = (accuracies{2}.accuracies);
values(values<58) = NaN;
scatter3(channel_sets(:, 1), channel_sets(:, 2), channel_sets(:, 3), values/10, values, 'filled');
xlabel('channel'); ylabel('channel'); zlabel('channel'); title('classification accuracy 3ch');
axis([0 16 0 16 0 16]);
caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
colorbar;
caxis([min(values) max(values)]);
view([4 25]); % [azimuth elevation]

%% Phi-star classification ACROSS FLIES 4D (4-channel) plot

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
load([results_directory results_filename]);

values = (accuracies{3}.accuracies);
channel_sets = double(accuracies{3}.channel_sets);

figure;
colormap('jet');
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices);
    values(values<58) = NaN;
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), values_sliced, values_sliced, '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    if fourth_d == min(channels_4d)
        colorbar;
    end
    axis([0 16 0 16 0 16]);
    caxis([58 max(accuracies{3}.accuracies)]);
    %caxis([min(accuracies{3}.accuracies) max(accuracies{3}.accuracies)]);
    %caxis([min(values) max(values)]);
    view([6 25]); % [azimuth elevation]
end

%% Phi-star classification ACROSS FLIES bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average classification accuracy across channel sets which include the channel

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
load([results_directory results_filename]);

figure;
colours = 'rgb';
for nChannels_counter = 1 : length(accuracies)
    nChannels = accuracies{nChannels_counter}.nChannels;
    channel_sets = double(accuracies{nChannels_counter}.channel_sets);
    values = accuracies{nChannels_counter}.accuracies;
    
    % Sum accuracies for each channel (sum across networks which contain the channel)
    values_collapsed = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
            end
        end
    end
    
    % Average accuracies
    values_collapsed = values_collapsed ./ set_counters;
    % Standard error
    values_collapsed_deviations = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
            end
        end
    end
    values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
    values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;
    
    % Plot
    %subplot(1, length(accuracies), nChannels_counter);
    plot((1:length(values_collapsed)), values_collapsed, colours(nChannels_counter)); hold on;
    %bar((1:length(values_collapsed)), values_collapsed);
    %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr);
    axis([0 16 50 75]);
    legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
    ylabel('%'); xlabel('channel'); title('Mean accuracy across sets containing channel X');
end

%% LOAD Phi-star values GLOBAL/NON-GLOBAL FLIES 2D, 3D, and 4D plots LOAD

% Global TPM (18s) was used in across fly analysis
% Non-global TPM (2.25s per trial) was used in within fly analysis
% All classification analyses were conducted at tau = 4ms

global_cov = 0; % 0: covariance over 2.25s; 1: covariance of 18s

if global_cov == 1
    % For across fly analysis we want one value per fly, so we use global covariance (first trial)
    disp('loading');
    data_directory = 'results/';
    data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar.mat'];
    load([data_directory data_filename]);
    disp('loaded');
else % global_cov == 0
    % For within fly analysis we use 16 values per fly (8 awake and 8 anest), so we used 2.25s TPMs
    disp('loading');
    data_directory = 'results/';
    data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar.mat'];
    load([data_directory data_filename]);
    disp('loaded');
end
tau = 1;

%% Phi-star values ACROSS FLIES 2D plot

% Build visual matrix
channel_sets = phis{1}.channel_sets;
value_map = zeros(max(channel_sets(:)));
values = permute(mean(mean(phis{1}.phi_stars(:, :, :, :, tau), 2), 3), [1 4 2 3]);
for value = 1 : size(values, 1)
    value_map(channel_sets(value, 1), channel_sets(value, 2)) = values(value, 1) - values(value, 2);
end
figure;
colormap('jet');
imagesc(value_map); cbar = colorbar; ylabel(cbar, '\Deltaphi3');
xlabel('channel'); ylabel('channel');
title('Awake - Anest');

%% Phi-star values ACROSS FLIES 3D plot

values = permute(mean(mean(phis{2}.phi_stars(:, :, :, :, tau), 2), 3), [1 4 2 3]);
channel_sets = phis{2}.channel_sets;
figure;
scatter3(channel_sets(:, 1), channel_sets(:, 2), channel_sets(:, 3), 75, values(:, 1) - values(:, 2), '.');
cbar = colorbar; ylabel(cbar, '\Deltaphi3');
xlabel('channel'); ylabel('channel'); zlabel('channel');
axis([0 16 0 16 0 16]);
view([4 25]); % [azimuth elevation]

%% Phi-star values 4D plot

values = permute(mean(mean(phis{3}.phi_stars(:, :, :, :, tau), 2), 3), [1 4 2 3]);
channel_sets = double(phis{3}.channel_sets);

% 4D to 2D PCA space projection
[coeff, score, latent] = pca(channel_sets);
figure;
colormap('jet');
scatter(score(:, 1), score(:, 2), [], (1:size(channel_sets, 1)), '.'); colorbar; % Map
figure;
colormap('jet');
plot_values = values(:, 1) - values(:, 2); % awake - anest
%plot_values = values(:, 1); % awake
%plot_values = values(:, 2); % anest
value_sizes = (rescale(plot_values, 1, 2000));
scatter(score(:, 1), score(:, 2), value_sizes, (1:size(channel_sets, 1)), '.'); colorbar;
set(gca,'Color','k'); % black background

% 4D+1 (4D space + phi measurement dimension) to 2D PCA projection
% Perhaps not so good because the spatial dimensions could be prioritised over the phi dimension
phi_dimension = (values(:, 1) - values(:, 2));
%channel_sets_scaled = rescale(channel_sets, min(phi_dimension), max(phi_dimension));
[coeff, score, latent, ~, explained] = pca(zscore([channel_sets_scaled phi_dimension]));
figure;
colormap('jet');
scatter(score(:, 1), score(:, 2), [], (1:size(channel_sets, 1)), '.'); colorbar;

% 4D to 1D PCA projection against phi-dimension
[coeff, score, latent] = pca(channel_sets);
phi_dimension = values(:, 1) - values(:, 2); % awake - anest
phi_dimension = values(:, 1);%values(:, 1); % awake
%phi_dimension = values(:, 2); % anest
figure;
colormap('jet');
scatter(score(:, 1), phi_dimension, [], (1:size(channel_sets, 1)), '.'); colorbar;
%set(gca,'Color','k'); % black background

% AWAKE - ANEST
figure;
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices, :);
    
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), 75, values_sliced(:, 1) - values_sliced(:, 2), '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    axis([0 16 0 16 0 16]);
    caxis([min(min(values(:, 1) - values(:, 2))) max(max(values(:, 1) - values(:, 2)))]);
    view([6 25]); % [azimuth elevation]
end

% AWAKE
figure;
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices, :);
    
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), 75, values_sliced(:, 1), '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    axis([0 16 0 16 0 16]);
    caxis([min(values(:, 2)) max(values(:, 1))]);
    view([6 25]); % [azimuth elevation]
end

% ANEST
figure;
channels_4d = channel_sets(:, 4);
for fourth_d = min(channels_4d) : max(channels_4d)
    subplot(3, 4, fourth_d + 1 - min(channels_4d));
    
    slice_logical_indices = channels_4d == fourth_d;
    
    channel_sets_sliced = channel_sets(slice_logical_indices, :);
    values_sliced = values(slice_logical_indices, :);
    
    scatter3(channel_sets_sliced(:, 1), channel_sets_sliced(:, 2), channel_sets_sliced(:, 3), 75, values_sliced(:, 2), '.');
    xlabel('channel'); ylabel('channel'); zlabel('channel'); title(fourth_d);
    axis([0 16 0 16 0 16]);
    caxis([min(values(:, 2)) max(values(:, 1))]);
    view([6 25]); % [azimuth elevation]
end

%% Phi-star values bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average phi value across channel sets which include the channel

figure;
colours = 'rgb';
for nChannels_counter = 1 : length(phis)
    nChannels = phis{nChannels_counter}.nChannels;
    channel_sets = double(phis{nChannels_counter}.channel_sets);
    values = permute(mean(mean(phis{nChannels_counter}.phi_stars(:, :, :, :, tau), 2), 3), [1 4 2 3]);
    values = values(:, 1) - values(:, 2);
    
    % Sum accuracies for each channel (sum across networks which contain the channel)
    values_collapsed = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
            end
        end
    end
    
    % Average phis
    values_collapsed = values_collapsed ./ set_counters;
    % Standard error
    values_collapsed_deviations = zeros(1, max(channel_sets(:)));
    set_counters = zeros(1, max(channel_sets(:)));
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                set_counters(channel) = set_counters(channel) + 1;
                values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
            end
        end
    end
    values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
    values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;
    
    % Plot
    %subplot(1, length(accuracies), nChannels_counter);
    plot((1:length(values_collapsed)), values_collapsed, colours(nChannels_counter)); hold on;
    %bar((1:length(values_collapsed)), values_collapsed);
    %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr, colours(nChannels_counter)); hold on;
    axis([0 16 0 0.03]);
    legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
    ylabel('\Phi'); xlabel('channel'); title('Mean phi3 across sets containing channel X');
end

%% Phi-star values WITHIN FLIES bar graph collapse

% Bar graph:
%   x-axis: channel sets which contain this channel
%   y-axis: average phi value across channel sets which include the channel

figure;
colours = 'rgb';
for fly = 1 : 13
    subplot(4, 4, fly);
    for nChannels_counter = 1 : length(phis)
        nChannels = phis{nChannels_counter}.nChannels;
        channel_sets = double(phis{nChannels_counter}.channel_sets);
        values = permute(mean(mean(phis{nChannels_counter}.phi_stars(:, :, fly, :, tau), 2), 3), [1 4 2 3]);
        values = values(:, 2);% - values(:, 2);
        
        % Sum accuracies for each channel (sum across networks which contain the channel)
        values_collapsed = zeros(1, max(channel_sets(:)));
        set_counters = zeros(1, max(channel_sets(:)));
        for channel = 1 : max(channel_sets(:))
            for channel_set = 1 : size(channel_sets, 1)
                if any(channel_sets(channel_set, :) == channel)
                    set_counters(channel) = set_counters(channel) + 1;
                    values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
                end
            end
        end
        
        % Average phis
        values_collapsed = values_collapsed ./ set_counters;
        % Standard error
        values_collapsed_deviations = zeros(1, max(channel_sets(:)));
        set_counters = zeros(1, max(channel_sets(:)));
        for channel = 1 : max(channel_sets(:))
            for channel_set = 1 : size(channel_sets, 1)
                if any(channel_sets(channel_set, :) == channel)
                    set_counters(channel) = set_counters(channel) + 1;
                    values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
                end
            end
        end
        values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
        values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;
        
        % Plot
        %subplot(1, length(accuracies), nChannels_counter);
        plot((1:length(values_collapsed)), values_collapsed, colours(nChannels_counter)); hold on;
        %bar((1:length(values_collapsed)), values_collapsed);
        %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr, colours(nChannels_counter)); hold on;
        axis([0 16 0 0.02]);
        legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
        ylabel('\Phi'); xlabel('channel'); title('Mean phi3 across sets containing channel X');
    end
end

%% Phi-star values ACROSS FLIES 4ch vs mean of channel means

% Get average 2-channel phi for each channel
nChannels_counter = 2;
nChannels = phis{nChannels_counter}.nChannels;
channel_sets = double(phis{nChannels_counter}.channel_sets);
values = permute(mean(mean(phis{nChannels_counter}.phi_stars(:, :, :, :, tau), 2), 3), [1 4 2 3]);
values = values(:, 1) - values(:, 2);

% Sum accuracies for each channel (sum across networks which contain the channel)
values_collapsed = zeros(1, max(channel_sets(:)));
set_counters = zeros(1, max(channel_sets(:)));
for channel = 1 : max(channel_sets(:))
    for channel_set = 1 : size(channel_sets, 1)
        if any(channel_sets(channel_set, :) == channel)
            set_counters(channel) = set_counters(channel) + 1;
            values_collapsed(channel) = values_collapsed(channel) + values(channel_set);
        end
    end
end

% Average phis
values_collapsed = values_collapsed ./ set_counters;
% Standard error
values_collapsed_deviations = zeros(1, max(channel_sets(:)));
set_counters = zeros(1, max(channel_sets(:)));
for channel = 1 : max(channel_sets(:))
    for channel_set = 1 : size(channel_sets, 1)
        if any(channel_sets(channel_set, :) == channel)
            set_counters(channel) = set_counters(channel) + 1;
            values_collapsed_deviations(channel) = values_collapsed_deviations(channel) + (values(channel_set)-values_collapsed(channel))^2;
        end
    end
end
values_collapsed_std = sqrt(values_collapsed_deviations ./ (set_counters-1));
values_collapsed_stderr = values_collapsed_std ./ set_counters.^0.5;

figure;
colours = 'rb';
shapes = 'ox';
subplot_counter = 1;
for nChannels_counter = 2 : length(phis)
    subplot(1, 2, nChannels_counter - 1);
    for condition = 1 : 2
        nChannels = phis{nChannels_counter}.nChannels;
        channel_sets = double(phis{nChannels_counter}.channel_sets);
        values = mean(mean(phis{nChannels_counter}.phi_stars(:, 1, :, condition, tau), 2), 3);
        
        % Find mean of mean 2-channel phis for each network
        values_2ch_mean = zeros(size(values));
        for channel_set = 1 : size(channel_sets, 1)
            channels = channel_sets(channel_set, :);
            values_2ch_mean(channel_set) = mean(values_collapsed(channels));
        end
        
        % Plot
        scatter(values_2ch_mean, values, [], colours(condition), shapes(condition)); hold on;
        disp([num2str(nChannels_counter+1) 'Ch, condition ' num2str(condition)]);
        [r, p] = corr(values_2ch_mean, values)
    end
    title([num2str(nChannels) ' channels']);
    xlabel('mean 2-ch phi3 across channels in set');
    ylabel('phi3 value');
    legend('awake', 'anest', 'Location', 'northeastoutside');
end
%% Description

%{

Plot classification performance of each concept order

All classification results are from using 4-channel phi (i.e. N = 1365)
All classification results are for classifying across flies (using
trial-averaged values)

%}

%% Setup

const_type = 'unpart'; % 'part'; 'both'
class_type = 'across'; % 'within'; 'across'

%% Load accuracies

results_dir = 'results/';
results_file = ['4ch_phi3Concept_' const_type '_svm_' class_type '.mat'];
acc = load([results_dir results_file]);

% Accuracies of individual concepts
c_accuracies = permute(mean(acc.net_accuracies(:, :, 1:15), 2), [1 3 2]);

% Accuracy of IIS
comp_accuracies = permute(mean(acc.net_accuracies(:, :, 16), 2), [1 3 2]);

% Accuracy of SII
p_accuracies = permute(mean(acc.net_accuracies(:, :, 17), 2), [1 3 2]);

%% Processing
% Average across all concepts within the same order

max_order = 4;

o_accuracies = zeros(size(c_accuracies, 1), max_order);

concept_ind = 1;
for order_c = 1 : max_order % for each order
    nConcepts = nchoosek(max_order, order_c);
    concept_indices = (concept_ind : concept_ind+nConcepts-1);
    o_accuracies(:, order_c) = mean(c_accuracies(:, concept_indices), 2); % mean across concepts
    %o_accuracies(:, order_c) = max(c_accuracies(:, concept_indices), [], 2); % max across concepts
    concept_ind = concept_ind + nConcepts;
end

%% Plot raincloud (for all; single concept, full constellation, big phi)

% Join together with big phi accuracies
acc_all = [o_accuracies comp_accuracies p_accuracies];
%acc_all = [normrnd(1, 0.5, [1000 1]) normrnd(2, 0.5, [1000 1]) normrnd(5, 0.5, [1000, 1])];

% Get correlation results
%tmp = load('../svm_classification_composition/results/4ch_medianSplit_correlation_svm_across.mat');
%r_acc = tmp.accuracy;
%acc_all = [r_acc' acc_all];

% How to set multiple baselines - https://stackoverflow.com/questions/44195924/mutiple-area-with-different-baseline-in-matlab

figure;
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
set(gcf, 'RendererMode', 'manual');
set(gcf, 'Renderer', 'painters');

colours = cbrewer('qual', 'Dark2', size(acc_all, 2)); % https://au.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/34087/versions/2/screenshot.jpg

base_offset = 37; % 17;
baselines = (0:-base_offset:-(size(acc_all, 2)-1)*base_offset); % baselines of clouds
rain_spread = 10;
cloud_rain_dist = 1;
rain_offset = (rain_spread/2) + cloud_rain_dist; % middle of rain
rain_scatter = (rand(size(acc_all, 1), 1) - 0.5) * rain_spread;

handles = cell(size(acc_all, 2), 1);

for feature_type = 1 : size(acc_all, 2)
    ax(feature_type) = axes; % Not possible to set multiple baselines within same axis;
    set(gca, 'XAxisLocation', 'top');
    
    if feature_type > 1
        ax(feature_type).YLim = ax(1).YLim;
        ax(feature_type).XLim = ax(1).XLim;
        hold on; % this locks the axes limits
    end
    
    % Create raincloud
    handles{feature_type} = raincloud_plot(acc_all(:, feature_type), 'box_on', 1, 'color', colours(feature_type, :), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    
    % Shift baseline (cloud)
    handles{feature_type}{1}.BaseValue = baselines(feature_type); % move baseline
    handles{feature_type}{1}.YData = handles{feature_type}{1}.YData + baselines(feature_type); % move cloud
    handles{feature_type}{1}.ShowBaseLine = 'off'; % hide baseline
    handles{feature_type}{1}.EdgeAlpha = 0; % hide cloud outline
    
    % Shift baseline (rain)
    handles{feature_type}{2}.YData = rain_scatter;
    handles{feature_type}{2}.YData = handles{feature_type}{2}.YData + baselines(feature_type) - rain_offset; % move rain
    handles{feature_type}{2}.SizeData = 2; % raindrop size
    
    % Shift boxplots
    handles{feature_type}{3}.Position([2 4]) = [-rain_spread/2 rain_spread]; % set box width to match rain spread
    handles{feature_type}{4}.YData = [-rain_spread/2 rain_spread/2]; % set median line to match rain spread
    handles{feature_type}{3}.Position(2) = handles{feature_type}{3}.Position(2) + baselines(feature_type) - rain_offset; % move box
    handles{feature_type}{4}.YData = handles{feature_type}{4}.YData + baselines(feature_type) - rain_offset; % move median line
    handles{feature_type}{5}.YData = handles{feature_type}{5}.YData + baselines(feature_type) - rain_offset; % move top whisker
    handles{feature_type}{6}.YData = handles{feature_type}{6}.YData + baselines(feature_type) - rain_offset; % move bot whisker
    
    % Hide all axes except the first
    if feature_type > 1
        axis off; % hide all axes except the first
    end
    
    view([-90 90]);
    
end

% Set same axis limits for all axes
limits = cell2mat(get(ax,'YLim')); % get both axes limits
set(ax,'YLim',[min(limits(:)) max(limits(:))]); % set the same values for both axes
linkaxes(ax);

% Link view angle, axes limits
%linkprop(ax, {'CameraPosition','CameraUpVector','YLim','XLim'});

%view([-90 90]); % view([90 -90]); % Swap x and y axes

% Plot chance line
plot([0.5 0.5], [-200 200], 'k--');

% Formatting
ylim(ax(1), [min(handles{end}{2}.YData)-2 max(handles{1}{1}.YData)+1]);
xlim(ax(1), [min(acc_all(:))-0.01 max(acc_all(:))+0.01]);
title(ax(1), ['classification using single concept (' const_type ')']);
xlabel(ax(1), [class_type '-fly accuracy']);
ylabel(ax(1), 'order of \phi');
set(ax(1), 'YTick', fliplr(baselines), 'YTickLabel', fliplr({'r', '1', '2', '3', '4', 'all', '\Phi'}));

% Draw separating lines
for base = max_order+1 : length(baselines)
    ypos = baselines(base)+(base_offset/2);
    plot([0 1], [ypos ypos], 'k');
end

%% Print

if strcmp(class_type, 'within')
    figure_name = 'figures/fig4a_raw'; % 4a = within-fly; 4b = across-flies
elseif strcmp(class_type, 'across')
    figure_name = 'figures/fig4b_raw';
else
    figure_name = 'figures/most_recent';
end

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG

%% Stats (LME)

channel_sets = nchoosek((1:15), 4);

% Concatenate accuracies into one matrix (sets x measures)
values = cat(2, c_accuracies, comp_accuracies, p_accuracies);

% label feature type (mechanism order, big-phi, composition)
% 1-4 = mechanism order
% 5 = composition; 6 = big-phi
type_labels = [1 1 1 1 2 2 2 2 2 2 3 3 3 3 4 5 6];

% Create table
tic;
table_array = zeros(numel(values), 5);
row_counter = 1;
for measure = 1 : size(values, 2)
    for network = 1 : size(values, 1)
        row = [network mean(channel_sets(network, :)) measure type_labels(measure) values(network, measure)];
        table_array(row_counter, :) = row;
        row_counter = row_counter + 1;
    end
end
toc

acc_table = array2table(table_array, 'VariableNames', {'set', 'set_center', 'feature', 'feature_type', 'accuracy'});
acc_table.set = categorical(acc_table.set);
acc_table.feature = categorical(acc_table.feature);
acc_table.feature_type = categorical(acc_table.feature_type);

%% LME with post-hocs

% LME
model_full_spec = 'accuracy ~ feature_type + (1|set)';
model_full = fitlme(acc_table, model_full_spec);

% Null
model_null_spec = 'accuracy ~ 1 + (1|set)';
model_null = fitlme(acc_table, model_null_spec);
compare(model_null, model_full)

% Post-hoc pairwise comparisons
model_spec = 'accuracy ~ feature_type + (1|set)';
comparisons = nchoosek((min(type_labels(:)):max(type_labels(:))), 2);
comparison_stats = cell(max(type_labels), max(type_labels));
for comparison = 1 : size(comparisons, 1)
    comparison_table = acc_table(ismember(double(acc_table.feature_type), comparisons(comparison, :)), :);
    comparison_table.feature_type = nominal(comparison_table.feature_type);
    disp('=================================================================');
    disp(comparisons(comparison, :));
    comparison_stats{comparisons(comparison, 1), comparisons(comparison, 2)} = fitlme(comparison_table, model_spec);
    disp('=================================================================');
end

%% Load SII values

nChannels = 4;

source_dir = '../phi_3/results/';
source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels' num2str(nChannels) '_globalTPM0.mat'];
tic
disp('loading');
tmp = load([source_dir source_file]);
disp('loaded');
toc

phis = tmp.phis{1};

%% Set-center regression (LME) for SII

channel_sets = nchoosek((1:15), nChannels);

% SII values
condition = 1;
values = permute(mean(phis.phis(:, :, :, condition), 2), [1 3 2 4]); % mean across trials

% Create table
tic;
table_array = zeros(numel(values), 5);
row_counter = 1;
for repeat = 1 : size(values, 2)
    for network = 1 : size(values, 1)
        row = [network repeat mean(channel_sets(network, :)) channel_set_distance(channel_sets(network, :)) values(network, repeat)];
        table_array(row_counter, :) = row;
        row_counter = row_counter + 1;
    end
end
toc

acc_table = array2table(table_array, 'VariableNames', {'set', 'repeat', 'set_center', 'dist', 'x'});
acc_table.set = categorical(acc_table.set);
acc_table.repeat = nominal(acc_table.repeat);

acc_table.x = zscore(acc_table.x);
%acc_table.set_center = zscore(acc_table.set_center);

%% LME - does include network location explain more variance?

tic;

% LME
model_full_spec = 'x ~ set_center*dist + (1|repeat)';
model_full = fitlme(acc_table, model_full_spec);

model_null_spec = cell(1, 1);
model_null_spec{1} = 'x ~ 1 + (1|repeat)';
model_null_spec{2} = 'x ~ set_center + dist + (1|repeat)';
model_null_spec{3} = 'x ~ set_center + set_center:dist + (1|repeat)';
model_null_spec{4} = 'x ~ dist + set_center:dist + (1|repeat)';
model_nulls = cell(size(model_null_spec));
for null_model = 1 : length(model_null_spec)
    model_nulls{null_model} = fitlme(acc_table, model_null_spec{null_model});
    compare(model_nulls{null_model}, model_full)
end

toc

%% Plot model
% Assuming x ~ set_center * dist + (1|repeat)

x = acc_table.set_center;
y = acc_table.dist;
zF = double(model_full.Coefficients(1, 2)) +...
    double(model_full.Coefficients(2, 2)).*x +...
    double(model_full.Coefficients(3, 2)).*y +...
    double(model_full.Coefficients(4, 2)).*x.*y;

figure;
scatter3(x, y, zF, [], zF, 'filled');

%% Plot SII actual values on center-distance map

addpath('../phi_3/');
networks = nchoosek((1:15), 4);
condition_labels = {'wake', 'anest'};

[map_interp, map, centers, distances] = center_distance_map(networks(acc_table.set, :), acc_table.x);

figure;
colormap inferno
imagesc(map_interp); colorbar;
set(gca, 'YDir', 'normal');
title(['SII map ' condition_labels{condition}]); ylabel('channel distance'); xlabel('set center');

%% Set-center regression (LME) for classification accuracy

% Accuracy values
measure = 17; % 16 = IIS; 17 = SII
channel_sets = nchoosek((1:15), 4);
values = mean(acc.net_accuracies(:, :, measure), 2);

% Create table
tic;
table_array = zeros(numel(values), 5);
row_counter = 1;
for network = 1 : size(values, 1)
    for repeat = 1 : size(values, 2)
        row = [network repeat mean(channel_sets(network, :)) channel_set_distance(channel_sets(network, :)) values(network)];
        table_array(row_counter, :) = row;
        row_counter = row_counter + 1;
    end
end
toc

acc_table = array2table(table_array, 'VariableNames', {'set', 'repeat', 'center', 'dist', 'x'});
acc_table.set = categorical(acc_table.set);

acc_table.x = zscore(acc_table.x);
%acc_table.set_center = zscore(acc_table.set_center);

%% LME - does include network location explain more variance?

% LME
model_full_spec = 'x ~ center + dist + (1|set)';
model_full = fitlme(acc_table, model_full_spec);

% Null
model_null_spec = cell(1, 1);
model_null_spec{1} = 'x ~ dist + (1|set)';
model_null_spec{2} = 'x ~ center + (1|set)';
model_nulls = cell(size(model_null_spec));
for null_model = 1 : length(model_null_spec)
    model_nulls{null_model} = fitlme(acc_table, model_null_spec{null_model});
    compare(model_nulls{null_model}, model_full)
end

%% Plot classification accuracies on center-distance map

addpath('../phi_3/');
networks = nchoosek((1:15), 4);
measure_labels{16} = 'IIS'; measure_labels{17} = 'SII';

[map_interp, map, centers, distances] = center_distance_map(networks(acc_table.set, :), acc_table.x);

figure;
colormap inferno
imagesc(map_interp); colorbar;
set(gca, 'YDir', 'normal');
title(['class. map ' measure_labels{measure} ' ' class_type]); ylabel('channel distance'); xlabel('set center');

%% Plot SII and classification accuracies on center-distance map, one figure

addpath('../phi_3/');
networks = nchoosek((1:15), 4);

figure;
colormap inferno

% Plot SII values
condition_labels = {'wake', 'anest'};
values = permute(mean(mean(phis.phis(:, :, :, :), 2), 3), [1 3 2 4]); % mean across trials, and flies
clim = [min(values(:)) max(values(:))];
positions = [1 4];
for condition = 1 : size(values, 4)
    subplot(2, 3, positions(condition));
    values_cond = values(:, :, :, condition);
    [map_interp, map, centers, distances, centers_axis, distances_axis] = center_distance_map(...
        repmat(networks, [size(values_cond, 2) 1]),... % values are (sets x repeats)
        values_cond(:));
    nan_locations = isnan(map_interp);
    imagesc(map_interp, 'AlphaData', ~nan_locations, clim); colorbar;
    set(gca, 'YDir', 'normal');
    set(gca, 'XTick', find(ismember(centers_axis, linspace(1, 15, 3))), 'XTickLabel', linspace(1, 15, 3));
    set(gca, 'YTick', linspace(1, length(distances_axis), 4), 'YTickLabel', distances_axis(linspace(1, length(distances_axis), 4)));
    set(gca, 'colorscale', 'log');
    title(['SII map ' condition_labels{condition}]);
    ylabel('channel distance'); xlabel('set center');
end

% Load classification values for IIS and SII
class_types = {'within', 'across'};
accuracies = []; % sets x IIS/SII x within/across
for class_type = 1 : length(class_types)
    
    % Load
    results_file = ['4ch_phi3Concept_' const_type '_svm_' class_types{class_type} '.mat'];
    acc = load([results_dir results_file]);
    
    % Accuracy of IIS
    type_accuracies = permute(mean(acc.net_accuracies(:, :, 16), 2), [1 3 2]); % average across repeats
    
    % Concatenate accuracy of SII
    type_accuracies = cat(2, type_accuracies,...
        permute(mean(acc.net_accuracies(:, :, 17), 2), [1 3 2])); % average across repeats
    
    accuracies = cat(3, accuracies, type_accuracies);
    
end

% Plot classification values
measure_labels = {'IIS', 'SII'};
%clim = [min(accuracies(:)) max(accuracies(:))];
positions = [2 3 5 6]; % within-IIS, within-SII, across-IIS, across-SII
pos_counter = 1;
for class_type = 1 : length(class_types)
    tmp = accuracies(:, :, class_type);
    clim = [min(tmp(:)) max(tmp(:))];
    for measure = 1 : size(accuracies, 2)
        subplot(2, 3, positions(pos_counter));
        [map_interp, map, centers, distances, centers_axis, distances_axis] = center_distance_map(...
            networks,...
            accuracies(:, measure, class_type));
        nan_locations = isnan(map_interp);
        imagesc(map_interp, 'AlphaData', ~nan_locations, clim); colorbar;
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', find(ismember(centers_axis, linspace(1, 15, 3))), 'XTickLabel', linspace(1, 15, 3));
        set(gca, 'YTick', linspace(1, length(distances_axis), 4), 'YTickLabel', distances_axis(linspace(1, length(distances_axis), 4)));
        title([measure_labels{measure} ' ' class_types{class_type}]);
        ylabel('channel distance'); xlabel('set center');
        pos_counter = pos_counter + 1;
    end
end

%% Print

figure_name = 'figures/fig5_raw';

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG
%% Description

%{

Plot classification performance of each concept order

All classification results are from using 4-channel phi (i.e. N = 1365)
All classification results are for classifying within flies (averaging
across flies)

%}

%% Setup

const_types = {'unpart', 'part', 'both'};
const_types = const_types(1);
class_type = 'within';

%% Load accuracies

c_accuracies = struct(); % accuracies for each individual concept

results_dir = 'results/';
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    results_file = ['4ch_phi3Concept_' const_type '_svm_' class_type '.mat'];
    acc_tmp = load([results_dir results_file]);
    
    % Take the values at a specific cost
    cost_param = 1;
    c_accuracies.(const_type) = mean(acc_tmp.cost_accuracies(:, :, :, cost_param), 2); % average across flies
    
    % Take the maximum across costs
    %c_accuracies.(const_type) = max(acc_tmp.cost_accuracies, [], 3);
    
end

%% Processing

% Average across all concepts within the same order

o_accuracies = c_accuracies; % accuracies averaged across concepts with the same order

max_order = acc_tmp.nChannels; % from the last loaded file, but all should have the same value

for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    
    % (networks x orders)
    o_accuracies.(const_type) = zeros(size(c_accuracies.(const_type), 1), max_order);
    
    concept_ind = 1; % 0 because we will add counters to it
    for order_c = 1 : max_order % for each order
        
        nConcepts = nchoosek(max_order, order_c);
        concept_indices = (concept_ind : concept_ind+nConcepts-1);
        
        o_accuracies.(const_type)(:, order_c) = mean(c_accuracies.(const_type)(:, concept_indices), 2);
        %o_accuracies.(const_type)(:, order_c) = max(c_accuracies.(const_type)(:, concept_indices), [], 2);
        
        concept_ind = concept_ind + nConcepts;
        
    end
    
end

% Big phi
p_accuracies = c_accuracies.unpart(:, end); % big-phi is same for all constellation types

%% Add big Phi results

% % Load classification results from using big phi (need to redo, using SVM?)
% results_dir = '../svm_classification_composition/results/';
% results_file = ['4ch_phiStarComposition_svm_' class_type '.mat'];
% 
% value_accuracies = struct();
% acc_tmp = load([results_dir results_file]);
% 
% % Take the values at a specific cost
% cost_param = 1;
% p_accuracies = acc_tmp.cost_accuracies(cost_param, :);
% 
% % Take the maximum across costs
% %p_accuracies = max(acc_tmp.cost_accuracies, [], 1);

%% Add full constellation results

% Load classification results from using full small phi constellation

results_dir = '../svm_classification_composition/results/';

comp_accuracies = struct(); % accuracies from using all concepts
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    
    results_file = ['4ch_phi3Composition_' const_type '_svm_' class_type '.mat'];
    acc_tmp = load([results_dir results_file]);
    
    % Take the values at a specific cost
    cost_param = 1;
    comp_accuracies.(const_type) = mean(acc_tmp.cost_accuracies(cost_param, :, :), 3); % average across flies
    
    % Take the maximum across costs
    %comp_accuracies.(const_type) = max(acc_tmp.cost_accuracies, [], 1);
end

%% Plot (line with errorbars)

const_offsets = [-0.2 0 0.2];

figure;

for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    
    errorbar((1:max_order) + const_offsets(const_type_c), mean(o_accuracies.(const_type), 1), std(o_accuracies.(const_type), [], 1)); hold on;
end

set(gca, 'XTick', (1:max_order));
xlabel('concept order');
ylabel('accuracy');
legend(const_types, 'Location', 'northwest');
title('classification using single concept');

%% Plot boxplot (just for single concepts)

const_offsets = [-0.2 0 0.2];

% Each row is for each const_type, each column is for each concept order
xlocations = repmat(const_offsets', [1 max_order]) + repmat((1:max_order), [length(const_types) 1]);
xlocations = xlocations';

% Join all values into single matrix, one box per column
o_acc_all = repmat(zeros(size(o_accuracies.(const_types{1}))), [1 length(const_types)]);
cols = (1:size(o_accuracies.(const_types{1}), 2));
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    o_acc_all(:, cols) = o_accuracies.(const_type);
    cols = cols + (length(cols));
end

figure;
boxplot(o_acc_all, 'Whisker', 5, 'Positions', xlocations(:));

set(gca, 'XTick', (1:max_order), 'XTickLabel', (1:max_order));
xlabel('concept order');
ylabel('accuracy');
title(['classification using single concept' newline '(left to right: unpart, part, both)']);

%% Plot raincloud

% Convert networks to data series, const_type to N repeated measures groups,
% orders to M measurements, giving M x N cell array

rain_accuracies = cell(max_order, length(const_types));
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    for order_c = 1 : max_order
        rain_accuracies{order_c, const_type_c} = o_accuracies.(const_type)(:, order_c);
    end
end

colours = [1 0 0; 0 1 0; 0 0 1];

figure;
h = rm_raincloud(rain_accuracies, colours);
set(gca, 'YLim', [-0.3 1.6]);
title(['classification using single concept']);
for i = 1  : size(h.s, 1)
    for j = 1 : size(h.s, 2)
        h.s{i, j}.SizeData = 1;
    end
end


%% Plot boxplot (for all; single concept, full constellation, big phi)

const_offsets = [-0.2 0 0.2];

figure;
%subplot(5, 1, (1:4));

% Single concept results
% Each row is for each const_type, each column is for each concept order
xlocations = repmat(const_offsets', [1 max_order]) + repmat((1:max_order), [length(const_types) 1]);
xlocations = xlocations';
% Join all values into single matrix, one box per column
o_acc_all = repmat(zeros(size(o_accuracies.(const_types{1}))), [1 length(const_types)]);
cols = (1:size(o_accuracies.(const_types{1}), 2));
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    o_acc_all(:, cols) = o_accuracies.(const_type);
    cols = cols + (length(cols));
end
boxplot(o_acc_all, 'Whisker', 5, 'Positions', xlocations(:));

hold on;

% Full constellation results
xlocations = xlocations(end, :) + 1; % Add new plots to the right
% Join all values into single matrix, one box per column
comp_acc_all = zeros(size(comp_accuracies.(const_types{1}), 2), length(const_types));
for const_type_c = 1 : length(const_types)
    const_type = const_types{const_type_c};
    comp_acc_all(:, const_type_c) = comp_accuracies.(const_type);
    cols = cols + (length(cols));
end
boxplot(comp_acc_all, 'Whisker', 5, 'Positions', xlocations);   

% Big phi results
xlocations = mean(xlocations(end, :)) + 1; % Add new plots to the right
boxplot(p_accuracies, 'Whisker', 5, 'Positions', xlocations);

% Plot chance line
plot([0 max(xlocations)+1], [0.5 0.5], 'k--');

% Formatting
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'XTick', (1:max(xlocations)), 'XTickLabel', {'1', '2', '3', '4', 'all', '\Phi'});
xlabel('concept order');
ylabel('accuracy');
title(['classification using single concept' newline '(left to right: unpart, part, both)']);
xlim([0.5 max(xlocations)+0.5]); ylim([0.3 0.9]);

%% Plot boxplot (for all; single concept, full constellation, big phi)
% For only one const_type

const_type = const_types{1};

const_offsets = [-0.2 0 0.2];

figure;
%subplot(5, 1, (1:4));

% Single concept results
o_acc_all = o_accuracies.(const_type);

% Full constellation results
comp_acc_all = comp_accuracies.(const_type);

% Plot with big phi results
boxplot([o_acc_all comp_acc_all' p_accuracies'], 'Whisker', 5);

hold on;

% Plot chance line
plot([0 max(xlocations)+1], [0.5 0.5], 'k--');

% Formatting
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'XTick', (1:max(xlocations)), 'XTickLabel', {'1', '2', '3', '4', 'all', '\Phi'});
xlabel('concept order');
ylabel('accuracy');
title(['classification using single concept' newline '(unpart)']);
xlim([0.5 max(xlocations)+0.5]); ylim([0.3 0.9]);

%% Plot raincloud (for all; single concept, full constellation, big phi)
% For only one const_type
% using rm_raincloud.m

const_type = const_types{1};

const_offsets = [-0.2 0 0.2];

% Single concept results
o_acc_all = o_accuracies.(const_type);

% Full constellation results
comp_acc_all = comp_accuracies.(const_type);

% Join together with big phi accuracies
acc_all = [o_acc_all comp_acc_all' p_accuracies'];

% Place accuracies into cell array
rain_accuracies = cell(size(acc_all, 2), 1);
for c = 1 : size(acc_all, 2)
    rain_accuracies{c} = acc_all(:, c);
end

colours = [0 0 0];

figure;
h = rm_raincloud(rain_accuracies, colours);
set(gca, 'YLim', [-0.3 1.6]);
title(['classification using single concept']);
for i = 1  : size(h.s, 1)
    for j = 1 : size(h.s, 2)
        h.s{i, j}.SizeData = 2;
    end
end

%% Plot raincloud (for all; single concept, full constellation, big phi)
% For only one const_type

const_type = const_types{1};

% Single concept results
o_acc_all = o_accuracies.(const_type);

% Max instead of mean across concepts within order
% o_acc_all(:, 1) = max(c_accuracies.unpart(:, (1:4)), [], 2);
% o_acc_all(:, 2) = max(c_accuracies.unpart(:, (5:10)), [], 2);
% o_acc_all(:, 3) = max(c_accuracies.unpart(:, (11:14)), [], 2);

% Full constellation results
comp_acc_all = comp_accuracies.(const_type);

% Join together with big phi accuracies
acc_all = [o_acc_all comp_acc_all' p_accuracies];

% Get correlation results
%tmp = load('../svm_classification_composition/results/4ch_medianSplit_correlation_svm_across.mat');
%r_acc = tmp.accuracy;
%acc_all = [r_acc' acc_all];

% Filter by divergence
% divs = permute(mean(tpm_divergence, 2), [1 3 4 2]);
% thresh = repmat(prctile(divs, 10, 1), [size(divs, 1) ones(length(size(divs))-1, 1)']);
% acc_all(divs(:, 11, 1) > thresh(:, 11, 1), :) = NaN;

% How to set multiple baselines - https://stackoverflow.com/questions/44195924/mutiple-area-with-different-baseline-in-matlab

figure;
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
set(gcf, 'RendererMode', 'manual');
set(gcf, 'Renderer', 'painters');

colours = cbrewer('qual', 'Dark2', size(acc_all, 2)); % https://au.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/34087/versions/2/screenshot.jpg

base_offset = 17;
baselines = (0:-base_offset:-(size(acc_all, 2)-1)*base_offset); % baselines of clouds
rain_spread = 5;
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
title(ax(1), ['classification using single concept (unpart)']);
xlabel(ax(1), 'accuracy');
ylabel(ax(1), 'order of \phi');
set(ax(1), 'YTick', fliplr(baselines), 'YTickLabel', fliplr({'r', '1', '2', '3', '4', 'all', '\Phi'}));

% Draw separating lines
for base = max_order+1 : length(baselines)
    ypos = baselines(base)+(base_offset/2);
    plot([0 1], [ypos ypos], 'k');
end

%% Print

figure_name = '../phi_3/figures/fig4_raw';

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG

%% Stats (LME)

channel_sets = nchoosek((1:15), 4);

% Concatenate accuracies into one matrix (sets x measures)
values = c_accuracies.unpart;
values = cat(2, permute(values, [1 3 2]), comp_accuracies.unpart');

% label feature type (mechanism order, big-phi, composition)
% 1-4 = mechanism order
% 5 = composition; 6 = big-phi
type_labels = [1 1 1 1 2 2 2 2 2 2 3 3 3 3 4 6 5];

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

%% LME - does include network location explain more variance?

% LME
model_full_spec = 'accuracy ~ feature_type + set_center + (1|set)';
model_full = fitlme(acc_table, model_full_spec);

% Null
model_null_spec = 'accuracy ~ feature_type + (1|set)';
model_null = fitlme(acc_table, model_null_spec);
compare(model_null, model_full)

%% Set-center LME

feature_type = categorical(4); % 5 = composition; 6 = big-phi
feature_table = acc_table(acc_table.feature_type == feature_type, :);

% LME
model_full_spec = 'accuracy ~ set_center + (1|set)';
model_full = fitlme(feature_table, model_full_spec);

% Null
model_null_spec = 'accuracy ~ 1 + (1|set)';
model_null = fitlme(feature_table, model_null_spec);
compare(model_null, model_full)

[r, p] = corr(feature_table.set_center, feature_table.accuracy)

%% Set-center, set-distance plot

addpath('../phi_3/');

[foo, foo2] = center_distance_map(nchoosek((1:15), 4), acc_all(:, 1)); % average across mechanisms with the same order
[foo, foo2] = center_distance_map(nchoosek((1:15), 4), c_accuracies.unpart(:, 3)); % individual mechanisms
figure; pcolor(foo); colorbar; caxis([min(acc_all(:)) max(acc_all(:))]); colormap inferno

%% Channel specificity

channel_sets = nchoosek((1:15), 4);

channel_acc = zeros(15, 1);

for channel = 1 : length(channel_acc)
    relevant_sets = any(channel_sets == channel, 2);
    channel_acc(channel) = mean(acc_all(relevant_sets, 5)); % average across mechanisms
    channel_acc(channel) = mean(c_accuracies.unpart(relevant_sets, 16)); % individual mechanisms
end

figure; plot(channel_acc);
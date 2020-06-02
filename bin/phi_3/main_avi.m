%% DESCRIPTION

%{

Plots composition Hasse graph (phi-3 is z-axis)

Average across all flies, for non-global TPM

See figures/videos from http://www.eneuro.org/content/4/5/ENEURO.0085-17.2017

%}

%% Setup

tau = 4;

%% Setup

output_file = 'animations/composition_all';

marker_size = 500;
nChannels = 4;

%% Load

load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials_nComponents4_phithree_nChannels4_globalTPM0.mat');

tau_levels = find(phis{1}.taus == tau);

%% Get state-weighted compositions for all parameters

composition_phis = phis{1}.big_mips(:, :, :, :, :, :, :, tau_levels);
%nonzero = composition_phis(composition_phis~=0);
%composition_phis = log(composition_phis+min(nonzero));
%composition_phis = log(composition_phis+1);

% Weight by state occurences (multiply phi by number of times the state occurred)
for partitioned = 1 : 2
    for concept = 1 : 15
        composition_phis(:, partitioned, concept, :, :, :, :) = ...
            permute(composition_phis(:, partitioned, concept, :, :, :, :), [1 4 5 6 7 2 3]) .* ...
            double(phis{1}.state_counters(:, :, :, :, :, tau_levels));
    end
end

% Sum across states
composition_phis = permute(sum(composition_phis, 1), [2 3 4 5 6 7 1]);

% Divide by total number of states (for weighted average)
% Assumes equal number of samples for all parameters
composition_phis = composition_phis ./ sum(phis{1}.state_counters(:, 1, 1, 1, 1, tau_levels));

% Average across trials
composition_phis = mean(composition_phis, 4);

% concepts x channel-sets x trials x flies x conditions x part-types
dims = size(composition_phis);
compositions = zeros([dims(2:end) 3]);

% Unpartitioned
composition_unpart = permute(composition_phis(1, :, :, :, :, :), [2 3 4 5 6 7 1]);
compositions(:, :, :, :, :, 1) = composition_unpart;
% Partitioned
composition_part = permute(composition_phis(2, :, :, :, :, :), [2 3 4 5 6 7 1]);
compositions(:, :, :, :, :, 2) = composition_part;
% Unpartitioned - partitioned
composition_diff = permute(composition_phis(1, :, :, :, :, :) - composition_phis(2, :, :, :, :, :), [2 3 4 5 6 7 1]);
compositions(:, :, :, :, :, 3) = composition_diff;

% % Log transform (transform each composition type separately)
% for part_type = 1 : 2 % (full minus disconnected is already normal)
%     compositions(:, :, :, :, :, part_type) = log(compositions(:, :, :, :, :, part_type));
% end
% % Log transform for unpartitioned minus partitioned
% nonzero = compositions(:, :, :, :, :, 3); nonzero = nonzero(nonzero~=0);
% compositions(:, :, :, :, :, 3) = log(compositions(:, :, :, :, :, 3)+(2*abs(min(nonzero))));

% Average across concepts within concept-order
dims = size(compositions);
orders = zeros([nChannels dims(2:end)]);
order_concepts = {(1:4), (5:10), (11:14), (15)};
for c_order = 1 : length(order_concepts)
    orders(c_order, :, :, :, :, :) = mean(compositions(order_concepts{c_order}, :, :, :, :, :), 1);
end

% Average orders across trials, channel-sets (flies x orders x conditions x part-types)
orders(orders == 0) = nan; % To ignore 0 values for first-order
orders_setMean = nanmean(mean(orders, 3), 2);
%orders_setMean = mean(log(mean(orders, 3)), 2);
orders_setMean = permute(orders_setMean, [4 1 5 6 2 3]);
%orders_setMean = log(orders_setMean);

% Average orders across trials, flies (channel-sets x orders x conditions x part-types)
orders_flyMean = mean(mean(orders, 3), 4);
%orders_flyMean = mean(log(mean(orders, 3)), 4);
orders_flyMean = permute(orders_flyMean, [2 1 5 6 3 4]);
%orders_flyMean = log(orders_flyMean);

% Average phi across trials, channel-sets (flies x conditions)
phi_setMean = mean(mean(phis{1}.phis(:, :, :, :, tau_levels), 2), 1);
% nonzeros = phis{1}.phis(phis{1}.phis ~= 0);
% phi_setMean = mean(mean(log(phis{1}.phis+min(nonzeros)), 2), 1);
% %phi_setMean = log(mean(mean(phis{1}.phis, 2), 1));
% phi_setMean = mean(log(mean(phis{1}.phis, 2)), 1);
phi_setMean = permute(phi_setMean, [3 4 2 1 5]);
phi_setMean = reshape(phis{1}.phis(:, :, :, :, tau_levels), [8*13 2]);

% Average phi across trials, flies (channel-sets x conditions)
phi_flyMean = mean(mean(phis{1}.phis(:, :, :, :, tau_levels), 2), 3);
% nonzeros = phis{1}.phis(phis{1}.phis ~= 0);
% phi_flyMean = mean(mean(log(phis{1}.phis+min(nonzeros)), 2), 3);
% %phi_flyMean = log(mean(mean(phis{1}.phis, 2), 3));
% phi_flyMean = mean(log(mean(phis{1}.phis, 2)), 3);
phi_flyMean = permute(phi_flyMean, [1 4 2 3 5]);

% t-score across flies, between conditions (channel-sets x orders x part-types)
composition_anest_tscores = zeros([size(orders, 2) size(compositions, 1) size(orders, 6)]);
trial_avg = permute(mean(compositions, 3), [1 2 4 5 6 3]);
for concept = 1 : size(trial_avg, 1)
    for net = 1 : size(trial_avg, 2)
        for part_type = 1 : size(trial_avg, 5)
            w = squeeze(trial_avg(concept, net, :, 1, part_type));
            a = squeeze(trial_avg(concept, net, :, 2, part_type));
            [h, p, ci, stats] = ttest(w, a);
            composition_anest_tscores(net, concept, part_type) = stats.tstat;
        end
    end
end

% Average concept tscores within concept-order
orders_anest_tscores = zeros(nChannels, size(composition_anest_tscores, 1), 3);
for c_order = 1 : length(order_concepts)
    orders_anest_tscores(c_order, :, :) = squeeze(mean(composition_anest_tscores(:, order_concepts{c_order}, :), 2));
    orders_anest_tscores(c_order, :, :) = squeeze(max(composition_anest_tscores(:, order_concepts{c_order}, :), [], 2));
end
orders_anest_tscores = permute(orders_anest_tscores, [2 1 3]); % (channel-sets x orders x part-types)

%% Rainclouds for everything Figure 2020-04-23

figure;
condition_colours = {'r', 'b'};
ylims = [-3.9 -2.8; -3.9 -2.8; -4.8 -4.15]; % for log-transformed values
ylims = [0.02 0.07; 0.02 0.07; 0 0.007];
set(gcf, 'color', 'w');
condition_offsets = [-0.2 0.2];
subplots = zeros(3, 2);

% Phi raincloud
subplot(size(subplots, 1), size(subplots, 2), 1);
data = cell(1, 1);
data{1, 1} = phi_setMean(:, 1);
data{2, 1} = phi_setMean(:, 2);
h = rm_raincloud(data, [0 0 0]);
% Update rain
for c = 1 : numel(h.s)
    h.s{c}.SizeData = 5;
end
% Update mean-points
for c = 1 : numel(h.m)
    h.m(c).SizeData = 10;
end
%set(gca, 'XScale', 'log');
title('\Phi');
xlabel('\Phi', 'Rotation', 0);
set(gca, 'YTick', 1+condition_offsets, 'YTickLabel', {'wake', 'anest'});
set(gca, 'YTick', [h.m(2).YData h.m(1).YData], 'YTickLabel', fliplr({'wake', 'anest'}));
delete(h.l);
box on;

% Phi raincloud (wake/anest)
subplot(size(subplots, 1), size(subplots, 2), 2);
data = cell(1, 1);
data{1} = phi_setMean(:, 1) ./ phi_setMean(:, 2);
h = rm_raincloud(data, [0 0 0]);
% Update rain
for c = 1 : numel(h.s)
    h.s{c}.SizeData = 5;
end
% Update mean-points
for c = 1 : numel(h.m)
    h.m(c).SizeData = 10;
end
% Labels
title('\Phi wake / anest.');
xlabel('\Delta\Phi', 'Rotation', 0); % x-y are swapped in rm_raincloud
set(gca, 'YTick', []);
%set(gca, 'XScale', 'log');
%hold on;
%line([0 0], [-100 100], 'Color', 'k');
box on;
hold on; line([1 1], get(gca, 'YLim'), 'LineWidth', 2, 'Color', 'k');

% Mechanism raindcloud
subplot(size(subplots, 1), size(subplots, 2), [3 4]); axis off;
axes_pos = get(gca, 'Position');
values = orders_setMean(:, :, :, 1); % 1=unpart; 2=part; 3=diff

base_offset = 100;
baselines_wake = (0:base_offset:(size(values, 2)-1)*base_offset);
baselines_anest = baselines_wake + (base_offset/4);
baselines = [baselines_wake; baselines_anest];
baselines = fliplr(baselines(:)');
%baselines = fliplr((0:base_offset:(numel(values)/size(values, 1)-1)*base_offset)); % baselines of clouds
rain_spread = 10;
cloud_rain_dist = 5;
rain_offset = (rain_spread/2) + cloud_rain_dist; % middle of rain
rain_scatter = (rand(size(values, 1), 1) - 0.5) * rain_spread;
colours = [1 0 0; 0 0 1];

cloud_counter = 1;
for c_order = 1 : size(values, 2)
    for condition = 1 : size(values, 3)
        ax(cloud_counter) = axes('Position', axes_pos);
        
        % Create raincloud
        handles{cloud_counter} = raincloud_plot(...
            values(:, c_order, condition),...
            'box_on', 0,...
            'color', colours(condition, :),...
            'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', 0,...
            'dot_dodge_amount', 0,...
            'box_col_match', 0);
        
        % Flatten cloud
        handles{cloud_counter}{1}.YData = handles{cloud_counter}{1}.YData * 0.25;
        
        % Shift baseline (cloud)
        handles{cloud_counter}{1}.BaseValue = baselines(cloud_counter); % move baseline
        handles{cloud_counter}{1}.YData = handles{cloud_counter}{1}.YData + baselines(cloud_counter); % move cloud
        handles{cloud_counter}{1}.ShowBaseLine = 'off'; % hide baseline
        handles{cloud_counter}{1}.EdgeAlpha = 0; % hide cloud outline
        
        % Shift baseline (rain)
        handles{cloud_counter}{2}.YData = rain_scatter;
        handles{cloud_counter}{2}.YData = handles{cloud_counter}{2}.YData + baselines(cloud_counter) - rain_offset; % move rain
        handles{cloud_counter}{2}.SizeData = 10; % raindrop size
        
        % Hide all axes except the first
        if cloud_counter > 1
            axis off;
        end
        
        view([-90 90]); % flip x-y axes
        
        %set(gca, 'XScale', 'log');
        
        set(gca,...
            'YTick', (baselines_wake+baselines_anest)/2,...
            'YTickLabel', (length(baselines_wake):-1:1));
        
        ylabel('mechanism size');
        xlabel('\phi');
        
        cloud_counter = cloud_counter + 1;
    end
end
% Set same axis limits for all axes (surounding space the width of raincloud)
% Assumes first cloud is the upper/leftmost
% [smallest baseline + space; highest point of highest cloud + space]
limits = [min(baselines)-cloud_rain_dist-(2*rain_spread)-cloud_rain_dist max(handles{1}{1}.YData)+rain_spread];
set(ax,'YLim',[min(limits(:)) max(limits(:))]); % set the same values for both axes
linkaxes(ax);

% Plot order rainclouds
part_type = 1; % only plot for UP
subplot(size(subplots, 1), size(subplots, 2), [5 6]);
% Calculate t-scores
% Convert to cell array (orders x conditions), each cell holds 1xflies
data = cell(size(orders_setMean, 2), 1);%size(orders, 3));
for c_order = 1 : size(data, 1)
    for condition = 1 : 1%size(data, 2)
        
        % Ratio wake/anest
        data{c_order, condition} = squeeze(orders_setMean(:, c_order, condition, part_type)) ./ squeeze(orders_setMean(:, c_order, 2, part_type)); % wake/anest
        
        %(for log transformed values)
        %data{c_order, condition} = squeeze(orders_flyMean(:, c_order, condition, part_type)) - squeeze(orders_flyMean(:, c_order, 2, part_type)); % wake-anest
        
        % Average tscores across all mechanisms with the same order
        %data{c_order, condition} = squeeze(orders_anest_tscores(:, c_order, part_type));
        
        % All tscores for all mechanisms
        %tscores = composition_anest_tscores(:, order_concepts{c_order}, part_type);
        %data{c_order, condition} = tscores(:);
    end
end
% Average concept tscores within concept-order
orders_anest_tscores = zeros(nChannels, size(composition_anest_tscores, 1), 3);
for c_order = 1 : length(order_concepts)
    orders_anest_tscores(c_order, :, :) = squeeze(mean(composition_anest_tscores(:, order_concepts{c_order}, :), 2));
end
orders_anest_tscores = permute(orders_anest_tscores, [2 1 3]); % (channel-sets x orders x part-types)
h = rm_raincloud(data, [0 0 0]);
%h = rm_raincloud(data, [1 0 0; 0 0 1]);
% Update rain
for c = 1 : numel(h.s)
    h.s{c}.SizeData = 5;
end
% Update lines
for c = 1 : numel(h.l)
    h.l(c).LineWidth = 1;
end
% Update mean-points
for c = 1 : numel(h.m)
    h.m(c).SizeData = 10;
end
% Labels
title('full (F) wake / anest.');
xlabel('\Delta\phi', 'Rotation', 0); ylabel('mechanism size'); % x-y are swapped in rm_raincloud
%xlim([0 0.45]);
%xlim([1 1.63]); ylim([-7 63]);
box on;
%hold on;
%line([0 0], [-100 800], 'Color', 'k');

%% Print

figure_name = '../figures/fig4_raw';

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG

%% Check bimodal distribution for 1-channel mechanisms

channels = (1:15);
networks = nchoosek(channels, 4);

order = 1;
condition = 1;

tmp = orders(order, :, :, :, condition, 1);
phi_lim = [min(tmp(:)) max(tmp(:))];

figure;
for fly = 1 : size(orders, 4)
    
    subplot(4, 4, fly);
    
    values = permute(orders(order, :, :, fly, condition, 1),...
        [2 1 3 4 5 6]);
    
    for channel_c = 1 : length(channels)
        channel = channels(channel_c);
        
        % Identify networks with the channel
        includes_channel = any(networks == channel, 2);
        
        % Add channel to data
        data{channel_c} = values(includes_channel);
        
    end
    
    % Plot channel data
    h = rm_raincloud(data, [0 0 0]);
    
    % Update rain
    for c = 1 : numel(h.s)
        h.s{c}.SizeData = 1;
    end
    
    % Update mean-points
    for c = 1 : numel(h.m)
        h.m(c).SizeData = 1;
    end
    
    % Update lines
    for c = 1 : numel(h.l)
        h.l(c).LineWidth = 1;
    end
    
    title(['fly' num2str(fly)]);
    xlabel('phi')
    ylabel('channel');
    xlim(phi_lim);
    
end

%% Build table for stats

% composition_unpart | composition_part | composition_diff
values = composition_unpart;

channel_sets = double(phis{1}.channel_sets);

% Determine concept order of each mechanism
order_labels = zeros(size(values, 1), 1);
counter = 1;
for m_order = 1 : size(channel_sets, 2)
    n_mechanisms = nchoosek(size(channel_sets, 2), m_order);
    order_labels(counter : counter + n_mechanisms - 1) = m_order;
    counter = counter + n_mechanisms;
end

% Create table (~20 seconds)
tic;
table_array = zeros(numel(values) / size(values, 3), length(size(values))+1); % exclude trials as a dimension
row_counter = 1;
for condition = 1 : size(values, 5)
    for fly = 1 : size(values, 4)
        for network = 1 : size(values, 2)
            for mechanism = 1 : size(values, 1)
                value = mean(values(mechanism, network, :, fly, condition), 3); % mean across trials
                row = [order_labels(mechanism) network mean(channel_sets(network, :)) fly condition value];
                table_array(row_counter, :) = row;
                row_counter = row_counter + 1;
            end
        end
    end
end
toc

comp_table = array2table(table_array, 'VariableNames', {'order', 'set', 'set_center', 'fly', 'condition', 'phi'});

comp_table_continuous = comp_table;

comp_table.fly = categorical(comp_table.fly);
comp_table.order = categorical(comp_table.order);
comp_table.condition = categorical(comp_table.condition);

%% LME
% Omnibus tests

% Full models
model_spec = 'phi ~ condition + order + condition:order + (1|fly) + (1|fly:set)';
model_full = fitlme(comp_table, model_spec);

% Null models (~15 seconds)
tic;
model_null_spec = cell(3, 1);
model_null_spec{1} = 'phi ~ condition + order + (1|fly) + (1|fly:set)';
model_null_spec{2} = 'phi ~ condition + condition:order + (1|fly) + (1|fly:set)';
model_null_spec{3} = 'phi ~ order + condition:order + (1|fly) + (1|fly:set)';
model_nulls = cell(size(model_null_spec));
for null_model = 1 : length(model_null_spec)
    model_nulls{null_model} = fitlme(comp_table, model_null_spec{null_model});
    compare(model_nulls{null_model}, model_full)
end
toc

%% Post-hoc tests among mechanism orders (awake)

wake_table = comp_table(comp_table.condition == categorical(1), :);

% 1-2 1-3 1-4 2-3 2-4 3-4

model_spec = 'phi ~ order + (1|fly) + (1|fly:set)';

comparisons = nchoosek((1:4), 2);
for comparison = 1 : size(comparisons, 1)
    comparison_table = wake_table(ismember(double(wake_table.order), comparisons(comparison, :)), :);
    comparison_table.order = nominal(comparison_table.order);
    disp('=================================================================');
    disp(comparisons(comparison, :));
    fitlme(comparison_table, model_spec)
    disp('=================================================================');
end

%% LME (absolute/relative change)

a = comp_table_continuous.phi(comp_table_continuous.condition == 1);
b = comp_table_continuous.phi(comp_table_continuous.condition == 2);
diff_table = comp_table_continuous(comp_table_continuous.condition == 1, :);
diff_table.phi = a ./ b; % a./b for original values, or a-b for log-transformed

diff_table.fly = categorical(diff_table.fly);
diff_table.set = categorical(diff_table.set);
diff_table.order = categorical(diff_table.order);

diff_table = diff_table(~isnan(diff_table.phi) & ~isinf(diff_table.phi), :);

% Full models
model_spec = 'phi ~ order + (1|fly) + (1|fly:set)';
model_full = fitlme(diff_table, model_spec);

% Null models
tic;
model_null_spec = cell(1, 1);
model_null_spec{1} = 'phi ~ 1 + (1|fly) + (1|fly:set)';
model_nulls = cell(size(model_null_spec));
for null_model = 1 : length(model_null_spec)
    model_nulls{null_model} = fitlme(diff_table, model_null_spec{null_model});
    compare(model_nulls{null_model}, model_full)
end
toc

%% LME for set-center (absolute/relative change)

a = comp_table_continuous.phi(comp_table_continuous.condition == 1);
b = comp_table_continuous.phi(comp_table_continuous.condition == 2);
diff_table = comp_table_continuous(comp_table_continuous.condition == 1, :);
diff_table.phi = a ./ b;

diff_table.fly = categorical(diff_table.fly);
diff_table.set = categorical(diff_table.set);
diff_table.order = categorical(diff_table.order);

diff_table = diff_table(~isnan(diff_table.phi) & ~isinf(diff_table.phi), :);

diff_table = diff_table(diff_table.order == categorical(4), :);

% Full models
model_spec = 'phi ~ set_center + (1|fly) + (1|fly:set)';
model_full = fitlme(diff_table, model_spec);

% Null models
tic;
model_null_spec = cell(1, 1);
model_null_spec{1} = 'phi ~ 1 + (1|fly) + (1|fly:set)';
model_nulls = cell(size(model_null_spec));
for null_model = 1 : length(model_null_spec)
    model_nulls{null_model} = fitlme(diff_table, model_null_spec{null_model});
    compare(model_nulls{null_model}, model_full)
end
toc

%% Post-hoc tests among mechanism orders (change in phi)

% 1-2 1-3 1-4 2-3 2-4 3-4

model_spec = 'phi ~ order + (1|fly) + (1|fly:set)';

comparisons = nchoosek((1:4), 2);
for comparison = 1 : size(comparisons, 1)
    comparison_table = diff_table(ismember(double(diff_table.order), comparisons(comparison, :)), :);
    comparison_table.order = nominal(comparison_table.order);
    disp('=================================================================');
    disp(comparisons(comparison, :));
    fitlme(comparison_table, model_spec)
    disp('=================================================================');
end

%% Big phi LME

% composition_unpart | composition_part | composition_diff
values = phis{1}.phis(:, :, :, :, tau_levels);

channel_sets = double(phis{1}.channel_sets);

% Create table
tic;
table_array = zeros(numel(values) / size(values, 2), length(size(values))+2); % exclude trials as a dimension
row_counter = 1;
for condition = 1 : size(values, 4)
    for fly = 1 : size(values, 3)
        for network = 1 : size(values, 1)
            value = mean(values(network, :, fly, condition), 2);
            row = [network mean(channel_sets(network, :)) channel_set_distance(channel_sets(network, :)) fly condition value];
            table_array(row_counter, :) = row;
            row_counter = row_counter + 1;
        end
    end
end
toc

phi_table = array2table(table_array, 'VariableNames', {'set', 'center', 'dist', 'fly', 'condition', 'phi'});

phi_table.fly = categorical(phi_table.fly);
phi_table.condition = categorical(phi_table.condition);

phi_table.phi = zscore(phi_table.phi);
%phi_table.center = zscore(phi_table.center);
%phi_table.dist = zscore(phi_table.dist);

% Omnibus tests

% Full models
model_spec = 'phi ~ condition + (1|fly) + (1|fly:set)';
model_full = fitlme(phi_table, model_spec);

% Null models (~15 seconds)
tic;
model_null_spec = cell(1, 1);
model_null_spec{1} = 'phi ~ 1 + (1|fly) + (1|fly:set)';
model_nulls = cell(size(model_null_spec));
for null_model = 1 : length(model_null_spec)
    model_nulls{null_model} = fitlme(phi_table, model_null_spec{null_model});
    compare(model_nulls{null_model}, model_full)
end
toc

%% Set center LME

% Full models
model_spec = 'phi ~ condition + center + dist + condition:center + condition:dist + (1|fly) + (1|fly:set)';
model_full = fitlme(phi_table, model_spec);

% Null models
tic;
model_null_spec = cell(1, 1);
model_null_spec{1} = 'phi ~ center + dist + condition:center + condition:dist + (1|fly) + (1|fly:set)';
model_null_spec{2} = 'phi ~ condition + dist + condition:center + condition:dist + (1|fly) + (1|fly:set)';
model_null_spec{3} = 'phi ~ condition + center + condition:center + condition:dist + (1|fly) + (1|fly:set)';
model_null_spec{4} = 'phi ~ condition + center + dist + condition:dist + (1|fly) + (1|fly:set)';
model_null_spec{5} = 'phi ~ condition + center + dist + condition:center + (1|fly) + (1|fly:set)';
model_nulls = cell(size(model_null_spec));
for null_model = 1 : length(model_null_spec)
    model_nulls{null_model} = fitlme(phi_table, model_null_spec{null_model});
    compare(model_nulls{null_model}, model_full)
end
toc
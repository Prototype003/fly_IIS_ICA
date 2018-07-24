%% Description

%{

Stats for figure 4 - comparing power, coherence, phi3, and phistar

Table variables are as follows:
value: classification accuracy (or some other value)
measure: power/coherence/phi3/phi*
nChannels: number of channels used in set
centre: set center
distance: set path distance
fly_id: fly number

Model for within classification
value ~ measure + nChannels + centre + distance + centre*distance + (1|fly_id)

Model for across classification
value ~ measure + nChannels + centre + distance + center*distance

%}

%% Setup

bin_location = '../';

addpath(bin_location); % We need the funtions to get set centres and path distances

freq_range_w = (1:42); %(1:83); % corresponding to ~5Hz and ~10Hz, check the 'frequencies' vector
freq_range_a = (1:329); %(1:329)=0-5Hz; There are more frequency bins for the single large trial

results_directory = [bin_location 'workspace_results/'];

%% Within model (average per fly

disp('within');

% Build table
disp('building table');

% Power
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification.mat';
load([results_directory results_filename]);
values = permute(mean(accuracies(freq_range_w, :, :), 1), [2 3 1]); % channels x flies
value = zeros(size(values, 2), 1);
measure = value + 5;
nChannels = value + 1;
fly_id = value;
row_counter = 1;
for fly = 1 : size(values, 2)
    value(row_counter) = mean(values(:, fly));
    fly_id(row_counter) = fly;
    row_counter = row_counter + 1;
end
table_within = table(value, measure, nChannels, fly_id);


% Coherency
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';
load([results_directory results_filename]);
values = permute(mean(accuracies(freq_range_w, :, :), 1), [2 3 1]);
value = zeros(size(values, 2), 1);
measure = value + 1;
nChannels = value + 2;
centre = value;
distance = value;
fly_id = value;
row_counter = 1;
for fly = 1 : size(values, 2)
        value(row_counter) = mean(values(:, fly));
        fly_id(row_counter) = fly;
        row_counter = row_counter + 1;
end
table_within = [table_within; table(value, measure, nChannels, fly_id)];

% Phi-three
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
load([results_directory results_filename]);
measure_counter = 6;
for nChannels_counter = 1 : length(accuracies)
    values = accuracies{nChannels_counter}.accuracies;
    value = zeros(size(values, 2), 1);
    measure = value + measure_counter;
    nChannels = value + accuracies{nChannels_counter}.nChannels;
    fly_id = value;
    row_counter = 1;
    for fly = 1 : size(values, 2)
            value(row_counter) = mean(values(:, fly));
            fly_id(row_counter) = fly;
            row_counter = row_counter + 1;
    end
    table_within = [table_within; table(value, measure, nChannels, fly_id)];
    measure_counter = measure_counter + 1;
end


% Phi-star
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
load([results_directory results_filename]);
measure_counter = 2;    
for nChannels_counter = 1 : length(accuracies)
    values = accuracies{nChannels_counter}.accuracies;
    value = zeros(size(values, 2), 1);
    measure = value + measure_counter;
    nChannels = value + accuracies{nChannels_counter}.nChannels;
    fly_id = value;
    row_counter = 1;
    for fly = 1 : size(values, 2)
            value(row_counter) = mean(values(:, fly));
            fly_id(row_counter) = fly;
            row_counter = row_counter + 1;
    end
    table_within = [table_within; table(value, measure, nChannels, fly_id)];
    measure_counter = measure_counter + 1;
end

% LME model
disp('building model');
model_spec = 'value ~ measure + (1|fly_id)';

model_within = fitlme(table_within, model_spec);

% Null models
disp('building null models');

model_null_specs = {...
    'value ~ 1 + (1|fly_id)'...
    };

% Build null models
model_nulls_within = cell(length(model_null_specs), 1);
for null = 1 : length(model_null_specs)
    disp(['fitting null model: ' model_null_specs{null}]);
    model_nulls_within{null} = fitlme(table_within, model_null_specs{null});
    disp('null model built');
end

% Model comparisons
disp('comparing models')
model_comparisons_within = cell(length(model_nulls_within), 1);
for null = 1 : length(model_nulls_within)
    disp(['null model: ' model_null_specs{null}]);
    model_comparisons_within{null} = compare(model_nulls_within{null}, model_within);
    disp(model_comparisons_within{null});
    disp('');
end

% Plot residuals
foo = fitted(model_within);
roo = residuals(model_within);
figure; % All measures
gscatter(foo, roo, table_within.measure);
% figure; % Power
% gscatter(foo(1:195), roo(1:195), table_within.measure(1:195));
% figure; % Coherence
% gscatter(foo(196:1560), roo(196:1560), table_within.measure(196:1560));
% figure; % Phi-3
% gscatter(foo(1561:26585), roo(1561:26585), table_within.measure(1561:26585));
% figure; % Phi-*
% gscatter(foo(26586:end), roo(26586:end), table_within.measure(26586:end));

% % Scatter plot for main effects
% figure;
% scatter(table_across.centre, table_across.value);


%% Within model

disp('within');

% Build table
disp('building table');

% Power
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification.mat';
load([results_directory results_filename]);
values = permute(mean(accuracies(freq_range_w, :, :), 1), [2 3 1]); % channels x flies
value = zeros(numel(values), 1);
measure = value + 5;
nChannels = value + 1;
centre = value;
distance = value;
fly_id = value;
row_counter = 1;
for fly = 1 : size(values, 2)
    for channel = 1 : size(values, 1)
        value(row_counter) = values(channel, fly);
        centre(row_counter) = channel;
        fly_id(row_counter) = fly;
        row_counter = row_counter + 1;
    end
end
table_within = table(value, measure, nChannels, centre, distance, fly_id);


% Coherency
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';
load([results_directory results_filename]);
values = permute(mean(accuracies(freq_range_w, :, :), 1), [2 3 1]);
value = zeros(numel(values), 1);
measure = value + 1;
nChannels = value + 2;
centre = value;
distance = value;
centres = mean(networks, 2);
distances = channel_set_distances(networks);
fly_id = value;
row_counter = 1;
for fly = 1 : size(values, 2)
    for set_counter = 1 : size(values, 1)
        value(row_counter) = values(set_counter, fly);
        centre(row_counter) = centres(set_counter);
        distance(row_counter) = distances(set_counter);
        fly_id(row_counter) = fly;
        row_counter = row_counter + 1;
    end
end
table_within = [table_within; table(value, measure, nChannels, centre, distance, fly_id)];

% Phi-three
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
load([results_directory results_filename]);
measure_counter = 6;
for nChannels_counter = 1 : length(accuracies)
    values = accuracies{nChannels_counter}.accuracies;
    value = zeros(numel(values), 1);
    measure = value + measure_counter;
    nChannels = value + accuracies{nChannels_counter}.nChannels;
    centre = value;
    distance = value;
    centres = mean(accuracies{nChannels_counter}.channel_sets, 2);
    distances = channel_set_distances(accuracies{nChannels_counter}.channel_sets);
    fly_id = value;
    row_counter = 1;
    for fly = 1 : size(values, 2)
        for set_counter = 1 : size(values, 1)
            value(row_counter) = values(set_counter, fly);
            centre(row_counter) = centres(set_counter);
            distance(row_counter) = distances(set_counter);
            fly_id(row_counter) = fly;
            row_counter = row_counter + 1;
        end
    end
    table_within = [table_within; table(value, measure, nChannels, centre, distance, fly_id)];
    measure_counter = measure_counter + 1;
end


% Phi-star
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
load([results_directory results_filename]);
measure_counter = 2;
for nChannels_counter = 1 : length(accuracies)
    values = accuracies{nChannels_counter}.accuracies;
    value = zeros(numel(values), 1);
    measure = value + measure_counter;
    nChannels = value + accuracies{nChannels_counter}.nChannels;
    centre = value;
    distance = value;
    centres = mean(accuracies{nChannels_counter}.channel_sets, 2);
    distances = channel_set_distances(accuracies{nChannels_counter}.channel_sets);
    fly_id = value;
    row_counter = 1;
    for fly = 1 : size(values, 2)
        for set_counter = 1 : size(values, 1)
            value(row_counter) = values(set_counter, fly);
            centre(row_counter) = centres(set_counter);
            distance(row_counter) = distances(set_counter);
            fly_id(row_counter) = fly;
            row_counter = row_counter + 1;
        end
    end
    table_within = [table_within; table(value, measure, nChannels, centre, distance, fly_id)];
    measure_counter = measure_counter + 1;
end

% LME model
disp('building model');
model_spec = 'value ~ measure + centre + distance + centre:distance + (1|fly_id)';

model_within = fitlme(table_within, model_spec);

% Null models
disp('building null models');

model_null_specs = {...
    'value ~ nChannels + centre + distance + centre:distance + (1|fly_id)',... % measure null model
    'value ~ measure + centre + distance + centre:distance + (1|fly_id)',... % nChannels null model
    'value ~ measure + nChannels + distance + centre:distance + (1|fly_id)',... % centre null model
    'value ~ measure + nChannels + centre + centre:distance + (1|fly_id)',... % distance main effect null model
    'value ~ measure + nChannels + centre + distance + (1|fly_id)'... % centre-distance interaction effect null model
    };

% Build null models
model_nulls_within = cell(length(model_null_specs), 1);
for null = 1 : length(model_null_specs)
    disp(['fitting null model: ' model_null_specs{null}]);
    model_nulls_within{null} = fitlme(table_within, model_null_specs{null});
    disp('null model built');
end

% Model comparisons
disp('comparing models')
model_comparisons_within = cell(length(model_nulls_within), 1);
for null = 1 : length(model_nulls_within)
    disp(['null model: ' model_null_specs{null}]);
    model_comparisons_within{null} = compare(model_nulls_within{null}, model_within);
    disp(model_comparisons_within{null});
    disp('');
end

% Plot residuals
foo = fitted(model_within);
roo = residuals(model_within);
figure; % All measures
gscatter(foo, roo, table_within.measure);
% figure; % Power
% gscatter(foo(1:195), roo(1:195), table_within.measure(1:195));
% figure; % Coherence
% gscatter(foo(196:1560), roo(196:1560), table_within.measure(196:1560));
% figure; % Phi-3
% gscatter(foo(1561:26585), roo(1561:26585), table_within.measure(1561:26585));
% figure; % Phi-*
% gscatter(foo(26586:end), roo(26586:end), table_within.measure(26586:end));

% % Scatter plot for main effects
% figure;
% scatter(table_across.centre, table_across.value);

%% Post-hoc tests (within)
% We will use signrank tests, averaging across channels/sets to obtain a single value per fly
% Assumed measure order:
%   1 - coherence
%   2, 3, 4 - phi-star(2, 3, 4)
%   5 - power
%   6 7 8 - phi-three(2, 3, 4)

table_summary = grpstats(table_within, {'measure', 'fly_id'});

% Power vs Phi-3
disp('power vs phi-3 (2ch, 3ch, 4ch)');
a = table_summary.mean_value(table_summary.measure == 5);
b = table_summary.mean_value(table_summary.measure == 6); % 2ch
[p, h, stats] = signrank(a, b);
disp(['p=' num2str(p) ' signedrank=' num2str(stats.signedrank)]);
b = table_summary.mean_value(table_summary.measure == 7); % 3ch
[p, h, stats] = signrank(a, b);
disp(['p=' num2str(p) ' signedrank=' num2str(stats.signedrank)]);
b = table_summary.mean_value(table_summary.measure == 8); % 4ch
[p, h, stats] = signrank(a, b);
disp(['p=' num2str(p) ' signedrank=' num2str(stats.signedrank)]);
disp(' ');

% Power vs Phi-*
disp('power vs phi-* (2ch, 3ch, 4ch)');
a = table_across.value(table_across.measure == 5);
b = table_across.value(table_across.measure == 2); % 2ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 3); % 3ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 4); % 4ch
[p, h, stats] = ranksum(b, a);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
disp(' ');

% Coherence vs Phi-3
disp('coherence vs phi-3 (2ch, 3ch, 4ch)');
a = table_across.value(table_across.measure == 1);
b = table_across.value(table_across.measure == 6); % 2ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 7); % 3ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 8); % 4ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
disp(' ');

% Coherence vs Phi-*
disp('power vs phi-* (2ch, 3ch, 4ch)');
a = table_across.value(table_across.measure == 1);
b = table_across.value(table_across.measure == 2); % 2ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 3); % 3ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 4); % 4ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
disp(' ');


%% Across model

disp('across');

% Build table
disp('building table');

% Power
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across1.mat';
load([results_directory results_filename]);
values = permute(mean(accuracies(freq_range_a, :), 1), [2 1]); % channels
value = zeros(numel(values), 1);
measure = value + 5;
nChannels = value + 1;
centre = value;
distance = value;
row_counter = 1;
for channel = 1 : size(values, 1)
    value(row_counter) = values(channel);
    centre(row_counter) = channel;
    row_counter = row_counter + 1;
end
table_across = table(value, measure, nChannels, centre, distance);


% Coherency
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification_across1.mat';
load([results_directory results_filename]);
values = permute(mean(accuracies(freq_range_a, :), 1), [2 1]);
value = zeros(numel(values), 1);
measure = value + 1;
nChannels = value + 2;
centre = value;
distance = value;
centres = mean(networks, 2);
distances = channel_set_distances(networks);
row_counter = 1;
for set_counter = 1 : size(values, 1)
    value(row_counter) = values(set_counter);
    centre(row_counter) = centres(set_counter);
    distance(row_counter) = distances(set_counter);
    row_counter = row_counter + 1;
end
table_across = [table_across; table(value, measure, nChannels, centre, distance)];

% Phi-three
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
load([results_directory results_filename]);
measure_counter = 6;
for nChannels_counter = 1 : length(accuracies)
    values = accuracies{nChannels_counter}.accuracies;
    value = zeros(numel(values), 1);
    measure = value + measure_counter;
    nChannels = value + accuracies{nChannels_counter}.nChannels;
    centre = value;
    distance = value;
    centres = mean(accuracies{nChannels_counter}.channel_sets, 2);
    distances = channel_set_distances(accuracies{nChannels_counter}.channel_sets);
    row_counter = 1;
    for set_counter = 1 : size(values, 1)
        value(row_counter) = values(set_counter);
        centre(row_counter) = centres(set_counter);
        distance(row_counter) = distances(set_counter);
        row_counter = row_counter + 1;
    end
    table_across = [table_across; table(value, measure, nChannels, centre, distance)];
    measure_counter = measure_counter + 1;
end

% Phi-star
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
load([results_directory results_filename]);
measure_counter = 2;
for nChannels_counter = 1 : length(accuracies)
    values = accuracies{nChannels_counter}.accuracies;
    value = zeros(numel(values), 1);
    measure = value + measure_counter;
    nChannels = value + accuracies{nChannels_counter}.nChannels;
    centre = value;
    distance = value;
    centres = mean(accuracies{nChannels_counter}.channel_sets, 2);
    distances = channel_set_distances(accuracies{nChannels_counter}.channel_sets);
    row_counter = 1;
    for set_counter = 1 : size(values, 1)
        value(row_counter) = values(set_counter);
        centre(row_counter) = centres(set_counter);
        distance(row_counter) = distances(set_counter);
        row_counter = row_counter + 1;
    end
    table_across = [table_across; table(value, measure, nChannels, centre, distance)];
    measure_counter = measure_counter + 1;
end

% LME model
disp('building model');
model_spec = 'value ~ measure + centre + distance + centre:distance';

model_across = fitlme(table_across, model_spec);

% Null models
disp('building null models');

model_null_specs = {...
    'value ~ nChannels + centre + distance + centre:distance',... % measure null model
    'value ~ measure + centre + distance + centre:distance',... % nChannels null model
    'value ~ measure + nChannels + distance + centre:distance',... % centre null model
    'value ~ measure + nChannels + centre + centre:distance',... % distance main effect null model
    'value ~ measure + nChannels + centre + distance'... % centre-distance interaction effect null model
    };

% Build null models
model_nulls_across = cell(length(model_null_specs), 1);
for null = 1 : length(model_null_specs)
    disp(['fitting null model: ' model_null_specs{null}]);
    model_nulls_across{null} = fitlme(table_across, model_null_specs{null});
    disp('null model built');
end

% Model comparisons
disp('comparing models')
model_comparisons_within = cell(length(model_nulls_across), 1);
for null = 1 : length(model_nulls_across)
    disp(['null model: ' model_null_specs{null}]);
    model_comparisons_within{null} = compare(model_nulls_across{null}, model_across);
    disp(model_comparisons_within{null});
    disp('');
end

% Plot residuals
foo = fitted(model_across);
roo = residuals(model_across);
figure; % All measures
gscatter(foo, roo, table_across.measure);
% figure; % Power
% gscatter(foo(1:15), roo(1:15), table_across.measure(1:15));
% figure; % Coherence
% gscatter(foo(16:120), roo(16:120), table_across.measure(16:120));
% figure; % Phi-3
% gscatter(foo(121:2045), roo(121:2045), table_across.measure(121:2045));
% figure; % Phi-*
% gscatter(foo(2046:end), roo(2046:end), table_across.measure(2046:end));

% % Scatter plot for main effects
% figure;
% scatter(table_across.centre, table_across.value);

%% Post-hoc tests (across)
% We will use ranksum tests, treating the pool of channels/sets as independent
% Assumed measure order:
%   1 - coherence
%   2, 3, 4 - phi-star(2, 3, 4)
%   5 - power
%   6 7 8 - phi-three(2, 3, 4)

% Power vs Phi-3
disp('power vs phi-3 (2ch, 3ch, 4ch)');
a = table_across.value(table_across.measure == 5);
b = table_across.value(table_across.measure == 6); % 2ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 7); % 3ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 8); % 4ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
disp(' ');

% Power vs Phi-*
disp('power vs phi-* (2ch, 3ch, 4ch)');
a = table_across.value(table_across.measure == 5);
b = table_across.value(table_across.measure == 2); % 2ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 3); % 3ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 4); % 4ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
disp(' ');

% Coherence vs Phi-3
disp('coherence vs phi-3 (2ch, 3ch, 4ch)');
a = table_across.value(table_across.measure == 1);
b = table_across.value(table_across.measure == 6); % 2ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 7); % 3ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 8); % 4ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
disp(' ');

% Coherence vs Phi-*
disp('power vs phi-* (2ch, 3ch, 4ch)');
a = table_across.value(table_across.measure == 1);
b = table_across.value(table_across.measure == 2); % 2ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 3); % 3ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
b = table_across.value(table_across.measure == 4); % 4ch
[p, h, stats] = ranksum(a, b);
disp(['p=' num2str(p) ' zval=' num2str(stats.zval) ' ranksum=' num2str(stats.ranksum)]);
disp(' ');

% Phi-3 nChannels comparison
% disp('phi-3 2v3v4ch');


% Phi-* nChannels comparison



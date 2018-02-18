%% Description

%{
Calculations are done in main_phi3_conditionalDependence.py
This script plots/does stats
%}

%% Setup

setSize = 4; % the set size where we want to find the maximum difference
setSize = setSize-1; 

data_directory = 'analysis_results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree_tpmDependence.mat';

%% Load

load([data_directory data_file]);

%% Plot for largest tpm with largest difference

% Find maximum difference
[max_diff, max_diff_ind] = max(tpms{setSize}.diff(:));
[set, fly, condition, tau] = ind2sub(size(tpms{setSize}.diff), max_diff_ind);

% Find maximum probability (for colour axis)
sbs1 = tpms{setSize}.sbs1(:, :, set, fly, condition, tau);
sbs2 = tpms{setSize}.sbs2(:, :, set, fly, condition, tau);
max_prob = max([max(sbs1(:)) max(sbs2(:))]);

sbn = tpms{setSize}.sbn(:, :, set, fly, condition, tau);

figure;
clim_sbs = [0 max_prob];
clim_sbn = [0 max(sbn(:))];
subplot(1, 5, [1 2]);
imagesc(tpms{setSize}.sbs1(:, :, set, fly, condition, tau), clim_sbs);
xlabel('state'); ylabel('state');
xticks([]); yticks([]);
c = colorbar;
title(c, 'p');

subplot(1, 5, 3);
imagesc(tpms{setSize}.sbn(:, :, set, fly, condition, tau), clim_sbn);
xlabel('channel'); ylabel('state');
xticks([]); yticks([]);
c = colorbar;
title(c, 'p');

subplot(1, 5, [4 5]);
imagesc(tpms{setSize}.sbs2(:, :, set, fly, condition, tau), clim_sbs);
xlabel('state'); ylabel('state');
xticks([]); yticks([]);
c = colorbar;
title(c, 'p');

%% Stats comparing differences between set sizes
% Need to use LME because number of sets changes, and sets observations are
% not independent

metric = 'corr'; % 'diff' or 'corr'

% Storage table
table_length = 0;
for setSize_counter = 1 : length(tpms)
    table_length = table_length + numel(tpms{setSize_counter}.(metric));
end
table_raw = struct();
table_raw.setSize = zeros(table_length, 1);
table_raw.set = zeros(table_length, 1);
table_raw.fly = zeros(table_length, 1);
table_raw.condition = zeros(table_length, 1);
table_raw.tau = zeros(table_length, 1);
table_raw.(metric) = zeros(table_length, 1);

table_row = 1;
set_label = 1;
tau_labels = [4 8 16];
for setSize_counter = 1 : length(tpms)
    metric_mat = tpms{setSize_counter}.(metric);
    
    if strcmp(metric, 'diff')
        % sum of differences divided by the number of the probabilities
        metric_mat = metric_mat / (2^(setSize_counter+1))^2;
    end
    
    for set = 1 : size(metric_mat, 1)
        for fly = 1 : size(metric_mat, 2)
            for condition = 1 : size(metric_mat, 3)
                for tau = 1 : size(metric_mat, 4)
                    table_raw.setSize(table_row) = setSize_counter + 1;
                    table_raw.set(table_row) = set;
                    table_raw.fly(table_row) = fly;
                    table_raw.condition(table_row) = condition;
                    table_raw.tau(table_row) = tau_labels(tau);
                    table_raw.(metric)(table_row) = fisher_rz(metric_mat(set, fly, condition, tau));
                    table_row = table_row + 1;
                end
            end
        end
        set_label = set_label + 1;
    end

end

metric_table = struct2table(table_raw);

%% Full model

model_spec = [metric ' ~ setSize + condition + tau + (1|fly) + (1|fly:set)'];
disp(['fitting model: ' model_spec]);
model_full = fitlme(metric_table, model_spec);
disp('full model built');

%% Null models
model_null_specs = {...
    [metric ' ~ setSize + tau + (1|fly) + (1|fly:set)'],...
    [metric ' ~ condition + tau + (1|fly) + (1|fly:set)'],...
    [metric ' ~ condition + setSize + (1|fly) + (1|fly:set)']...
    };

model_nulls = cell(length(model_null_specs), 1);
for null = 1 : length(model_null_specs)
    disp(['fitting null model: ' model_null_specs{null}]);
    model_nulls{null} = fitlme(metric_table, model_null_specs{null});
    disp('null model built');
end

%% Model Comparison

compare(model_nulls{2}, model_full);

%%

metric = 'corr';
tau = 3;
condition_colours = {'r', 'b'};

figure;

for fly = 1 : size(tpms{setSize_counter}.(metric), 2)
    for condition = 1 : size(tpms{setSize_counter}.(metric), 3)
        means = zeros(3, 1);
        stds = zeros(3, 1);
        for setSize_counter = 1 : length(tpms)
            metric_mat = tpms{setSize_counter}.(metric);
            if strcmp(metric, 'diff')
                metric_mat = metric_mat / (2^(setSize_counter+1)^2);
            end
            means(setSize_counter) = mean(metric_mat(:, fly, condition, tau));
            stds(setSize_counter) = std(metric_mat(:, fly, condition, tau)) / sqrt(size(metric_mat, 1));
        end
        errorbar([2 3 4], means, stds, condition_colours{condition});
        hold on;
    end
end
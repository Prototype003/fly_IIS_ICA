%% Description

%{
Conduct MDS on correlation of correlations
%}

%% Settings

flies = (1:13);

%% Load

% load fly data
load('../workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat');

% median split
% fly_meds = repmat(median(fly_data, 1), [size(fly_data, 1) 1 1 1 1]);
% fly_data = fly_data > fly_meds;

%% Correlation between channels

% channels x channels x trials x conditions x flies
channel_corrs = zeros(size(fly_data, 2), size(fly_data, 2), size(fly_data, 3), size(fly_data, 5), length(flies));

for fly = flies
    for trial = 1 : size(fly_data, 3)
        for condition = 1 : size(fly_data, 5)
            channel_corrs(:, :, trial, condition, fly) = corr(fly_data(:, :, trial, fly, condition));
        end
    end
end

channel_corrs = channel_corrs .^ 2;

%% Upper triangle only

upper_sample = triu(ones(size(channel_corrs(:, :, 1, 1, 1))), 1);
upper_inds = find(upper_sample);

channel_corrs_vec = zeros(length(upper_inds), size(channel_corrs, 3), size(channel_corrs, 4), size(channel_corrs, 5));

for fly = flies
    for trial = 1 : size(channel_corrs, 3)
        for condition = 1 : size(channel_corrs, 4)
            tmp = channel_corrs(:, :, trial, condition, fly);
            channel_corrs_vec(:, trial, condition, fly) = tmp(upper_inds);
        end
    end
end

%% Correlation distance matrix

% Reshape to collapse conditions into trials

trial_patterns = permute(channel_corrs_vec, [4 1 2 3]); % flies x corrs x trials x conditions
trial_patterns = reshape(trial_patterns, [size(trial_patterns, 1) size(trial_patterns, 2) size(trial_patterns, 3)*size(trial_patterns, 4)]);
trial_patterns = permute(trial_patterns, [2 3 1]); % corrs x trials x flies

trial_corrs = zeros(size(trial_patterns, 2), size(trial_patterns, 2), size(trial_patterns, 3));

for fly = flies
    trial_corrs(:, :, fly) = corr(trial_patterns(:, :, fly));
end

corr_dists = 1 - trial_corrs;
corr_dists = round(corr_dists, 10);

%% Correlation of correlation patterns

pattern_corrs = zeros(size(channel_corrs));

for fly = flies
    for trial = 1 : size(channel_corrs, 3)
        for condition = 1 : size(channel_corrs, 4)
            pattern_corrs(:, :, trial, condition, fly) = corr(channel_corrs(:, :, trial, condition, fly));
        end
    end
end

%% Upper triangle only

upper_sample = triu(ones(size(pattern_corrs(:, :, 1, 1, 1))), 1);
upper_inds = find(upper_sample);

pattern_corrs_vec = zeros(length(upper_inds), size(pattern_corrs, 3), size(pattern_corrs, 4), size(pattern_corrs, 5));

for fly = flies
    for trial = 1 : size(pattern_corrs, 3)
        for condition = 1 : size(pattern_corrs, 4)
            tmp = pattern_corrs(:, :, trial, condition, fly);
            pattern_corrs_vec(:, trial, condition, fly) = tmp(upper_inds);
        end
    end
end

%% Correlation distance matrix

% Reshape to collapse conditions into trials

trial_patterns = permute(pattern_corrs_vec, [4 1 2 3]); % flies x corrs x trials x conditions
trial_patterns = reshape(trial_patterns, [size(trial_patterns, 1) size(trial_patterns, 2) size(trial_patterns, 3)*size(trial_patterns, 4)]);
trial_patterns = permute(trial_patterns, [2 3 1]); % corrs x trials x flies

trial_corrs = zeros(size(trial_patterns, 2), size(trial_patterns, 2), size(trial_patterns, 3));

for fly = flies
    trial_corrs(:, :, fly) = corr(trial_patterns(:, :, fly));
end

corr_dists = 1 - trial_corrs;
corr_dists = round(corr_dists, 10);

%% MDS

upper_sample = triu(ones(size(corr_dists(:, :, 1))));
upper_inds = find(upper_sample);

scaled = cell(length(flies));

for fly = flies
    tmp = corr_dists(:, :, fly);
    scaled{fly} = cmdscale(tmp);
end

%% Plot

for fly = flies
    figure;
    scatter(scaled{fly}(1:8, 1), scaled{fly}(1:8, 2), 'ro'); hold on;
    scatter(scaled{fly}(9:16, 1), scaled{fly}(9:16, 2), 'bx');
end

%% Mean across flies

scaled = cmdscale(mean(corr_dists, 3));
figure;
scatter(scaled(1:8, 1), scaled(1:8, 2), 'ro'); hold on;
scatter(scaled(9:16, 1), scaled(9:16, 2), 'bx');

%% MDS across all flies and trials

data_all = permute(fly_data, [1 2 3 5 4]); % samples x channels x trials x conditions x flies
data_all = reshape(data_all, [size(data_all, 1) size(data_all, 2) size(data_all, 3)*size(data_all, 4) size(data_all, 5)]); % collapse conditions
data_all = reshape(data_all, [size(data_all, 1) size(data_all, 2) size(data_all, 3)*size(data_all, 4)]); % collapse across flies

% Final dimensions should be samples x channels x observations
%   where observations is trials nested in conditions nested in flies

% Correlation among channels at each trial
channel_corrs = zeros(size(data_all, 2), size(data_all, 2), size(data_all, 3));
for trial = 1 : size(data_all, 3)
    channel_corrs(:, :, trial) = corr(data_all(:, :, trial));
end
channel_corrs = channel_corrs .^ 2;

% Correlation among patterns of channel correlations
pattern_corrs = zeros(size(channel_corrs));
for trial = 1 : size(channel_corrs, 3)
    pattern_corrs(:, :, trial) = corr(channel_corrs(:, :, trial));
end

% Correlation distances among trials
upper_sample = triu(ones(size(pattern_corrs(:, :, 1))), 1);
upper_inds = find(upper_sample);
pattern_corrs_vec = zeros(length(upper_inds), size(pattern_corrs, 3));
for trial = 1 : size(pattern_corrs, 3)
    tmp = channel_corrs(:, :, trial);
    pattern_corrs_vec(:, trial) = tmp(upper_inds);
end
corr_dists = 1 - corr(pattern_corrs_vec);

% MDS
scaled = cmdscale(round(corr_dists, 9));

% Plot
figure;
cond_colours = repmat((1:size(fly_data, 5)), [size(fly_data, 3), 1]);
for fly = 1 : size(fly_data, 4)
    from = (fly-1)*size(fly_data, 3)*size(fly_data, 5) + 1;
    to = fly*size(fly_data, 3)*size(fly_data, 5);
    scatter(scaled((from:to), 1), scaled((from:to), 2), 100, cond_colours(:), '.'); hold on;
    text(scaled((from:to), 1), scaled((from:to), 2), [' ' num2str(fly)], 'FontSize', 6);
end
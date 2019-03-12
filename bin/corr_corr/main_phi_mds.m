%% Description

%{
Conduct MDS on correlation of correlations
%}

%% Settings

flies = (1:13);
nChannels = 2;
tau = 1;

%% Load

addpath('../figure_code/');

phis = phi_load('phi_three', 0, '../');

max_channel = max(phis{nChannels-1}.channel_sets(:));
phi_data = phis{nChannels-1}.phis(:, :, :, :, tau);

%% Collapse conditions and flies

% phi_data (sets x trials x flies x conditions)

channel_phis_vec = permute(phi_data(:, :, :, :), [1 2 4 3]); % sets x trials x conditions x flies
channel_phis_vec = reshape(channel_phis_vec, [size(channel_phis_vec, 1) size(channel_phis_vec, 2)*size(channel_phis_vec, 3) size(channel_phis_vec, 4)]); % collapse conditions
channel_phis_vec = reshape(channel_phis_vec, [size(channel_phis_vec, 1) size(channel_phis_vec, 2)*size(channel_phis_vec, 3)]); % collapse across flies

%% Treat phi as the correlation

% correlation of phi patterns
% convert phi vectors into "correlation matrices"

channel_phis = zeros([repmat(max_channel, [1 nChannels]) size(channel_phis_vec, 2)]);
for trial = 1 : size(channel_phis_vec, 2)
    tmp_vec = channel_phis_vec(:, trial);
    tmp_mat = zeros(repmat(max_channel, [1 nChannels]));
    tmp_mat(tril(true(repmat(max_channel, [1 nChannels])), -1)) = tmp_vec;
    channel_phis(:, :, trial) = tmp_mat' + tmp_mat;
end

% Correlation among patterns of channel correlations
pattern_corrs = zeros(size(channel_phis));
for trial = 1 : size(channel_phis, 3)
    pattern_corrs(:, :, trial) = corr(channel_phis(:, :, trial));
end

% Correlation distances among trials
upper_sample = triu(ones(size(pattern_corrs(:, :, 1))), 1);
upper_inds = find(upper_sample);
pattern_corrs_vec = zeros(length(upper_inds), size(pattern_corrs, 3));
for trial = 1 : size(pattern_corrs, 3)
    tmp = pattern_corrs(:, :, trial);
    pattern_corrs_vec(:, trial) = tmp(upper_inds);
end

corr_dists = 1 - corr(pattern_corrs_vec);
corr_dists = 1 - corr(channel_phis_vec);

% MDS
scaled = cmdscale(round(corr_dists, 9));

% If there is an opportunity to get the matrix to work

% Plot
figure;
cond_colours = repmat((1:size(fly_data, 5)), [size(fly_data, 3), 1]);
for fly = 1 : size(fly_data, 4)
    from = (fly-1)*size(fly_data, 3)*size(fly_data, 5) + 1;
    to = fly*size(fly_data, 3)*size(fly_data, 5);
    scatter(scaled((from:to), 1), scaled((from:to), 2), 100, cond_colours(:), '.'); hold on;
    text(scaled((from:to), 1), scaled((from:to), 2), [' ' num2str(fly)], 'FontSize', 6);
end

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
    tmp = pattern_corrs(:, :, trial);
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
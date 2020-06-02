%% Description

%{

Plots power spectrum of independent components

%}

%% Settings

% Number of a priori independent components we want to get
nComponents = 4;

%% Setup

% Sampling rate (Hz)
sample_rate = 1000;

% Load ICs
load(['results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials'...
    '_nComponents' num2str(nComponents)...
    '_all.mat']);
load(['results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials'...
    '_nComponents' num2str(nComponents)...
    '.mat']);

%% Get overall weights for the transformation (RICA)

ica_models = ica_models(1:13); % incorrect length specified when this was created

weights = zeros([size(ica_models{1}.TransformWeights) length(ica_models)]);
for fly = 1 : length(ica_models)
    weights(:, :, fly) = pca_models.coeffs(:, :, fly) * ica_models{fly}.TransformWeights;
end

weights = abs(weights); % To see which are most important

% Plot
nFlies = 13;
figure;
for fly = 1 : nFlies
    subplot(4, 4, fly);
    imagesc(weights(:, :, fly)', [min(weights(:)) max(weights(:))]);
    cbar = colorbar; title(cbar, 'abs. weight');
    ylabel('IC'); xlabel('channel');
    title(['fly' num2str(fly)]);
end

% Order ICs based on channel location

% Find max weight across channels for each IC
[value, location] = max(weights, [], 1);
% Sort ICs based on channel location of max weight
[sorted, order] = sort(location, 2);
% Reorder weights
for fly = 1 : size(weights, 3)
    weights(:, :, fly) = weights(:, order(1, :, fly), fly);
end

% Plot
nFlies = 13;
figure;
for fly = 1 : nFlies
    subplot(4, 4, fly);
    imagesc(weights(:, :, fly)', [min(weights(:)) max(weights(:))]);
    cbar = colorbar; title(cbar, 'abs. weight');
    ylabel('IC'); xlabel('channel');
    title(['fly' num2str(fly)]);
end

%% Get overall weights for the transformation (fastica)

weights = abs(ica_models.coeffs); % To see which are most important

% Plot
nFlies = 13;
figure;
for fly = 1 : nFlies
    subplot(4, 4, fly);
    imagesc(weights(:, :, fly)', [min(weights(:)) max(weights(:))]);
    cbar = colorbar; title(cbar, 'abs. weight');
    ylabel('IC'); xlabel('channel');
    title(['fly' num2str(fly)]);
end

% Order ICs based on channel location

% Find max weight across channels for each IC
[value, location] = max(weights, [], 1);
% Sort ICs based on channel location of max weight
[sorted, order] = sort(location, 2);
% Reorder weights
for fly = 1 : size(weights, 3)
    weights(:, :, fly) = weights(:, order(1, :, fly), fly);
end

% Plot
nFlies = 13;
figure;
for fly = 1 : nFlies
    subplot(4, 4, fly);
    imagesc(weights(:, :, fly)', [min(weights(:)) max(weights(:))]);
    cbar = colorbar; title(cbar, 'abs. weight');
    ylabel('IC'); xlabel('channel');
    title(['fly' num2str(fly)]);
end

%% Compute power spectra

L = size(fly_data, 1); % Length of signal

Y = fft(fly_data(:, :, :, :, :), [], 1);

% Two-sided spectrum
P2 = abs(Y/L);

% One-sided spectrum
P1 = P2(1:L/2+1, :, :, :, :);
P1(2:end-1, :, :, :, :) = 2*P1(2:end-1, :, :, :, :);

f = sample_rate*(0:(L/2))/L;

%% Plot power spectra

fly = 3;
trials = (8);

cond_colours = {'r', 'b'};

figure;
for condition = 1 : size(P1, 5)
    data = P1(:, :, trials, fly, condition);
    data = mean(data, 3); % average across trials
    data = reshape(data, [size(data, 1) numel(data)/size(data, 1)]);
    plot(f, data, cond_colours{condition}); hold on;
end
title('Single-Sided Amplitude Spectrum');
xlabel('f (Hz)');
ylabel('|P1(f)|');

%% Plot original signals

fly = 2;
trials = (8);

cond_colours = {'r:', 'b:'};

figure;
for condition = 1 : size(fly_data, 5)
    plot((1:size(fly_data, 1)), mean(fly_data(:, :, trials, fly, condition), 3), cond_colours{condition}); hold on;
end
xlim([1 size(fly_data, 1)]);
title('IC time-course');
xlabel('t');
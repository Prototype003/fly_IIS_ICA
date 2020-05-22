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
load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICA_all.mat');

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

fly = 1;
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

fly = 1;
trials = (8);

cond_colours = {'r:', 'b:'};

figure;
for condition = 1 : size(fly_data, 5)
    plot((1:size(fly_data, 1)), mean(fly_data(:, :, trials, fly, condition), 3), cond_colours{condition}); hold on;
end
xlim([1 size(fly_data, 1)]);
title('IC time-course');
xlabel('t');
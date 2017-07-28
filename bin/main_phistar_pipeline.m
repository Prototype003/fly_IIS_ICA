%% DESCRIPTION

%{
This is for creating the figure showing the phi-star calculation pipeline
(it visualises sample data, and the associated covariance matrices)

%}

%% Load
data_directory = 'workspace_results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';
load([data_directory '/' data_file '.mat']);

%% Compute covariance matrices
fly = 5;
condition = 2;
trial = 4;
channels = (11:14);
tau = 4;
[cov_present_present, cov_present_past, cov_past_past] = Cov_comp_sample(fly_data(:, channels, trial, fly, condition)', tau);

%% Plot
figure;

% Raw
sub1 = subplot(2, 3, (1:3));
imagesc(fly_data((101:120), channels, trial, fly, condition)');
ylabel('channel');
cbar = colorbar;
xlabel(cbar, '\muV');
colormap(sub1, 'jet');

% Covariances
subplot(2, 3, 4);
imagesc(cov_present_present);
title('cov(t, t)');
subplot(2, 3, 5);
imagesc(cov_present_past);
title('cov(t, t-\tau)');
subplot(2, 3, 6);
imagesc(cov_past_past);
title('cov(t-\tau, t-\tau)');
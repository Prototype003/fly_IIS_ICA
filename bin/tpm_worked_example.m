%{
Figure creation for research proposal

Experimental setup > data > binarised data > transition probability matrix
> phi
%}

%{
dimensions = [3 10]; % channels x samples
sample_range = [-30 60];

sample_data = rand(dimensions);
sample_data = sample_range(1) + (sample_range(2) - sample_range(1)) .* sample_data;

middle = median(sample_data(:));
middle_mat = repmat(middle, dimensions);

binarised_data = sample_data > middle_mat;

figure;
colormap('jet');
subplot(4, 1, [1 2]);
imagesc(sample_data);

subplot(4, 1, [3 4]);
imagesc(binarised_data); colorbar;
%}

%{
% Load
%load('../../flies/fly_data/trials_anesthDescAdded11092014_bPlrRerefTyp2/Analyzed_WalkingDrosDror115610072014/trials.mat');
load('../../flies/fly_data_lineNoiseRemoved/Analyzed_WalkingDrosDror115610072014/trials.mat');

channels = (1:2);
samples = (81787:81787+19);

sample_data = [trials.LFP];
sample_data = sample_data(channels, samples);

% Each channel is binarised based on its median value
middle = median(sample_data, 2);
middle_mat = repmat(middle, [1 length(samples)]);
%middle_mat = repmat(median(sample_data(:)), [length(channels) length(samples)]);

binarised_data = sample_data > middle_mat;

figure;
sub1 = subplot(4, 1, [1 2]);
imagesc(sample_data); cbar = colorbar;
set(gca, 'YTick', [1 2 3], 'XTickLabel', '');
xlabel(cbar, '\muV');
colormap(sub1, 'jet');

sub2 = subplot(4, 1, [3 4]);
imagesc(binarised_data);
set(gca, 'YTick', [1 2 3]);
cbar = colorbar; caxis([0 1]);
set(cbar, 'YTick', [0.25 0.75], 'YTickLabel', {'off', 'on'});
xlabel('time sample'); ylabel('channel');
colormap(sub2, [0 0 0; 1 1 1]);
%}

%% Settings

fly = 1;
channels = [5 6];
trial = 4;
samples = (1:10); %(101:110);
condition = 1;
tau = 1;

%% Load
% data_directory = 'workspace_results/';
% data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';
% 
% disp("Loading fly data");
% loaded_data = load([data_directory data_file '.mat']);
% fly_data = loaded_data.fly_data; % Reminder: dimensions are: (samples x channels x trials x flies x conditions)
% disp("Fly data loaded")

%%

% Binarise fly_data
%fly_data_binarised = binarise_global_median(fly_data);
n_values = 2;

% Get relevant data
raw_data = fly_data(samples, channels, trial, fly, condition);

% Each channel is binarised based on its median value
middle = median(raw_data, 1);
middle_mat = repmat(middle, [size(raw_data, 1), 1]);
binarised_data = raw_data > middle_mat;
channel_data = binarised_data; %fly_data_binarised(samples, channels, trial, fly, condition);

% Actual TPM
tpm = build_tpm(channel_data, 1, n_values);

% Independent TPM
tpm_ind = build_tpm_independent(channel_data, 1, n_values)

%% Plot data
figure;
data_plot = subplot(2, 1, 1);
imagesc(raw_data'); cbar = colorbar;
set(gca, 'YTick', [1 2 3], 'XTickLabel', '');
xlabel(cbar, '\muV');
colormap(data_plot, 'jet');

binarised_plot = subplot(2, 1, 2);
imagesc(channel_data');
set(gca, 'YTick', [1 2]);
cbar = colorbar; caxis([0 1]);
set(cbar, 'YTick', [0.25 0.75], 'YTickLabel', {'off', 'on'});
xlabel('time sample'); ylabel('channel');
colormap(binarised_plot, [0 0 0; 1 1 1]);

%% Plot TPMs

figure;
imagesc(tpm); colorbar;

%% Function: build independent TPM from 2 channel data
function [tpm] = build_tpm_independent(channel_data, tau, n_values)
% For the 2 channel scenario
% Multiplies two (independent) TPMs together using Kronecker Tensor
% multiplication (kron())
% Note that the output is NOT in the required LOLI format for PyPhi

% Build TPM for single channels
a = build_tpm(channel_data(:, 1), tau, n_values);
b = build_tpm(channel_data(:, 2), tau, n_values);

% Multiply single-channel TPMs
tpm = kron(a, b);

end
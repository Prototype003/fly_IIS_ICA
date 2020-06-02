%% Description

%{

Converts independent components back into channels

%}

%% Global settings

% Number of a priori independent components we want to get
nComponents = 4;

%% Setup

addpath('FastICA_25/');

% Load data
file_prefix = ['results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials_nComponents' num2str(nComponents)];
data_file = file_prefix;
load(data_file); % ICs
data_file = [file_prefix '_all']; % ICA model information
load(data_file);

% Output filename
out_prefix = [file_prefix '_ic2channels'];

%% Use channel weight matrices to re-obtain original dimensionality

ic_data = fly_data;
nChannels = size(ica_models.coeffs, 1);

% Flatten trials/conditions dimensions into time-samples dimension
ic_data = permute(ic_data, [1 3 5 2 4]);
dims = size(ic_data);
ic_data = reshape(ic_data, [prod(dims(1:3)) dims(4:end)]);

dims = size(ic_data); dims(2) = nChannels;
data = zeros(dims);

for fly = 1 : size(fly_data, 4)
    data(:, :, fly) = ic_data(:, :, fly) * ica_models.mixers(:, :, fly);
end

%% Convert back to original dimensions

dims = size(fly_data);
dims(2) = nChannels;

% Separate conditions
data = reshape(data, [dims(1)*dims(3) dims(5) dims(2) dims(4)]);

% Separate trials
data = reshape(data, [dims(1) dims(3) dims(5) dims(2) dims(4)]);

data = permute(data, [1 4 2 5 3]);

fly_data = data; % rewrite IC data

%% Save

save(out_prefix, 'fly_data');
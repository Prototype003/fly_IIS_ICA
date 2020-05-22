%% Description

%{

Conducts Independent Component Analysis on fly LFPs

Outputs independent components to results/

%}

%% Global settings

% Number of a priori independent components we want to get
nComponents = 4;

%% Setup

addpath('FastICA_25/');

% Load data
data_file = '../../../fly_phi/bin/workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat';
load(data_file);

% Output filename
out_prefix = ['results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials_nComponents' num2str(nComponents)];

%% Reformat data to avoid constantly having to use nested loops

% field_names = {'LFP', 'trial', 'fly', 'condition'};
% fly_table = array2table(zeros(0, 4), 'VariableNames', field_names);
% for fly = 1 : size(fly_data, 4)
%     for condition = 1 : size(fly_data, 5)
%         for trial = 1 : size(fly_data, 3)
%             fly_table = [fly_table; cell2table({fly_data(:, :, trial, fly, condition), trial, fly, condition}, 'VariableNames', field_names)];
%         end
%     end
% end

%% Combine trials/conditions for ICA

% samples x trials x conditions x channels x flies
data = permute(fly_data, [1 3 5 2 4]);

dims = size(data);

data = reshape(fly_data, [prod(dims(1:3)) dims(4) dims(5)]);

%% Conduct ICA
tic;

pc_scores = zeros(size(data));
pc_coeffs = zeros(size(fly_data, 2), size(fly_data, 2), size(fly_data, 4));
pc_latents = zeros(size(fly_data, 2), size(fly_data, 4));

ica_models = cell(size(data, 2), 1);
ica_scores = zeros(size(data, 1), nComponents, size(data, 3));

for fly = 1 : size(data, 3)
    disp(fly);
    
    % Whiten the data using zscore and PCA
    [coeff, score, latent] = pca(zscore(data(:, :, fly)));
    pc_scores(:, :, fly) = score;
    pc_coeffs(:, :, fly) = coeff;
    pc_latents(:, fly) = latent;
    
    % ICA
    ica_models{fly} = rica(score, nComponents); % MATLAB's reconstructionICA
    
    % Get independent components
    ic = transform(ica_models{fly}, score);
    ica_scores(:, :, fly) = ic;
    
end

toc

%% Convert IC scores to original format
% (original format - can reuse older code)
% NOTE - original fly_data variable is overwritten with ICs

dims = size(fly_data); % the second dimension is changed by ICA
dims(2) = nComponents;

fly_data = zeros(dims);

for row = 1 : size(fly_table, 1)
    
    fly_data(:, :, fly_table.trial(row), fly_table.fly(row), fly_table.condition(row)) = fly_ica.ic_score{row};
    
end

%% Save

% TODO - update to have number of independent components in the filename

% Save everything
save([out_prefix '_all.mat'], 'fly_data', 'fly_table', 'fly_pca', 'fly_rica', 'fly_ica');

% Python can't load MATLAB tables, save only IC matrix
save([out_prefix '.mat'], 'fly_data');
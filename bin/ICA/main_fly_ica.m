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
out_prefix = ['results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICA_nComponents' num2str(nComponents)];

%% Reformat data to avoid constantly having to use nested loops

field_names = {'LFP', 'trial', 'fly', 'condition'};
fly_table = array2table(zeros(0, 4), 'VariableNames', field_names);
for fly = 1 : size(fly_data, 4)
    for condition = 1 : size(fly_data, 5)
        for trial = 1 : size(fly_data, 3)
            fly_table = [fly_table; cell2table({fly_data(:, :, trial, fly, condition), trial, fly, condition}, 'VariableNames', field_names)];
        end
    end
end

%% Conduct ICA
tic;

pca_models = cell(size(fly_table, 1), 3);
rica_models = cell(size(fly_table, 1), 1);
independent_components = cell(size(fly_table, 1), 1);

for row = 1 : size(fly_table)
    disp(row);
    
    lfp = fly_table.LFP{row};
    
    % Whiten the data using zscore and PCA
    [coeff, score, latent] = pca(zscore(lfp));
    pca_models(row, :) = {coeff, score, latent};
    
    % ICA
    rica_models{row} = rica(score, nComponents); % MATLAB's reconstructionICA
    
    % Get independent components
    ic = transform(rica_models{row}, score);
    independent_components{row} = ic;
    
end

% Convert to table
fly_pca = cell2table(pca_models, 'VariableNames', {'pc_coeff', 'pc_score', 'pc_latent'});
fly_rica = cell2table(rica_models, 'VariableNames', {'rica_model'});
fly_ica = cell2table(independent_components, 'VariableNames', {'ic_score'});

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

% Save everything
save([out_prefix '_all.mat'], 'fly_data', 'fly_table', 'fly_pca', 'fly_rica', 'fly_ica');

% Python can't load MATLAB tables, save only IC matrix
save([out_prefix '.mat'], 'fly_data');
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

%% Combine trials/conditions for ICA

% samples x trials x conditions x channels x flies
data = permute(fly_data, [1 3 5 2 4]);

dims = size(data);

data = reshape(data, [prod(dims(1:3)) dims(4) dims(5)]);

%% Conduct ICA
tic;

pc_scores = zeros(size(data));
pc_coeffs = zeros(size(fly_data, 2), size(fly_data, 2), size(fly_data, 4));
pc_latents = zeros(size(fly_data, 2), size(fly_data, 4));

ica_models = cell(size(data, 3), 1);
ica_scores = zeros(size(data, 1), nComponents, size(data, 3));

ica_mixers = zeros(nComponents, size(data, 2), size(data, 3));
ica_coeffs = zeros(size(data, 2), nComponents, size(data, 3));

for fly = 1 : size(data, 3)
    disp(fly);
    
%     % Whiten the data using zscore and PCA
%     [coeff, score, latent] = pca(zscore(data(:, :, fly)));
%     pc_scores(:, :, fly) = score;
%     pc_coeffs(:, :, fly) = coeff;
%     pc_latents(:, fly) = latent;
%     
%     % ICA
%     ica_models{fly} = rica(score, nComponents); % MATLAB's reconstructionICA
%     
%     % Get independent components
%     % NOTE - final weights are just coeff * ica_model.TransformWeights
%     ic = transform(ica_models{fly}, score);
%     ica_scores(:, :, fly) = ic;
    
    % ICA
    % NOTE - fastica uses dimensions (signals x samples)
    % NOTE - (AB)' = B'A'
    % NOTE - refer to line 481 of fastica
    %   icasig = W * mixedsig = (signals x samples)
    %   (samples x signals) = (W*mixedsig)' = (mixedsig'*W')
    %       where mixedsig as input to fastica is signals x samples)
    [ics, mix, demix] = fastica(data(:, :, fly)', 'numOfIC', nComponents);
    ica_scores(:, :, fly) = ics';
    ica_mixers(:, :, fly) = mix';
    ica_coeffs(:, :, fly) = demix';
    
end

pca_models = struct();
pca_models.scores = pc_scores;
pca_models.coeffs = pc_coeffs;
pca_models.latents = pc_latents;

ica_models = struct();
ica_models.scores = ica_scores;
ica_models.coeffs = ica_coeffs;
ica_models.mixers = ica_mixers;

toc

%% Convert IC scores to original format
% (original format - can reuse older code)
% NOTE - original fly_data variable is overwritten with ICs

new_dims = dims;
new_dims(4) = nComponents; % number of channels has changed

ica_scores = reshape(ica_scores, new_dims);
ica_scores = permute(ica_scores, [1 4 2 5 3]);

fly_data = ica_scores;

%% Save

% Save everything
save([out_prefix '_all.mat'], 'ica_models');

% Python can't load MATLAB tables, save only IC matrix
save([out_prefix '.mat'], 'fly_data');
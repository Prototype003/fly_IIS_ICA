%% Description

%{

Demonstration of ICA

%}

%% Random sources (undirectly observed, original signals)

a = 0.5*sind((1:1000))';
b = (rand(1000, 1)-0.5);

figure; scatterhist(a, b);
figure;
subplot(2, 1, 1); plot(a);
subplot(2, 1, 2); plot(b);

%% Mixed sources (observed signals)
% Observed signals are linearly mixed
% Note:
%   If we have more observed signals than sources, model is over-determined
%       (which is OK)
%   If we have fewer observed signals than sources, model is under-
%       determined (which is generally not OK, but there are methods for
%       dealing with this, e.g. MATLAB's RICA)

x = 1*a + 5*b;
y = 2*a + 3*b;

figure; scatterhist(x, y);
figure;
subplot(2, 1, 1); plot(x);
subplot(2, 1, 2); plot(y);

%% Sphering - PCA
% Remove correlations between the observed signals
% Note - sphering (PCA and whitening) is done by fastica, which also
%   re-centers the ICs at the end

data = [x y]; % samples x signals

[coeff, score, latent] = pca(data);

figure; scatterhist(score(:, 1), score(:, 2));
figure;
subplot(2, 1, 1); plot(score(:, 1));
subplot(2, 1, 2); plot(score(:, 2));

%% Sphering - scaling

data_white = zscore(score); % 'white' - the process is called 'whitening'

figure; scatterhist(data_white(:, 1), data_white(:, 2));

%% ICA (FastICA)
% Goal of ICA is to minimise Gaussianity of the sources
% Note - FastICA takes matrix dimensions (signals x samples), opposite to
%   the general dimensions used for MATLAB functions (samples x signals)

addpath('FastICA_25/');

data = [x y];

data_ic = fastica(data');
data_ic = data_ic';

figure; scatterhist(data_ic(:, 1), data_ic(:, 2));
figure;
subplot(2, 1, 1); plot(data_ic(:, 1));
subplot(2, 1, 2); plot(data_ic(:, 2));

%% Scatter plot ICs
% Rotate using view([angle 90]);

figure; scatter(data_ic(:, 1), data_ic(:, 2));

%% MATLAB's RICA
% Note - for some reason, MATLAB's RICA gives different result to FastICA

rica_model = rica(data, 2); % We want 2 ICs

data_rica = transform(rica_model, data_white);

figure; scatterhist(data_rica(:, 1), data_rica(:, 2));
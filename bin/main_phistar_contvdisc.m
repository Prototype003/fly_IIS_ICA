%% DESCRIPTION

%{
This is for comparing continuous phi-star with discrete phi-star (they should be correlated)

%}

%% SETUP

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

filter_percent = 5; % 0 to 100 %

condition_shapes{1} = 'o'; condition_shapes{2} = 'x';
tau_colours = {'r', 'g', 'b'};
tau_alphas = [1 0.8 0.5];

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    ];

%% LOAD

% disp('loading');
% 
% % Continuous
% load([data_directory data_filename '_phistar.mat']);
% continuous = phis;
% 
% % Median split
% load([data_directory data_filename '_medianSplit1_phistar.mat']);
% discrete = phis;
% 
% disp('loaded');

%% Correlate all data points at each nChannels

figure;
for nChannels_counter = 1 : length(continuous)
    subplot(1, length(continuous), nChannels_counter);
    
    % Plot all data for each tau
    for tau_counter = 1 : size(continuous{nChannels_counter}.phi_stars, 5)
        continuous_points = continuous{nChannels_counter}.phi_stars(:, :, :, :, tau_counter);
        discrete_points = discrete{nChannels_counter}.phi_stars(:, :, :, :, tau_counter);
        scatter(continuous_points(:), discrete_points(:), 100, ['.' tau_colours{tau_counter}], 'MarkerEdgeAlpha', tau_alphas(tau_counter));
        hold on;
    end
end


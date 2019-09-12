%% 

%{
Plots datapoints which are used for classification
%}

%% SETUP

class_type = 'across'; % 'across' or 'within'

nFeatures = 1;

tau = 1;

% bin directory location
bin_location = '../';

addpath('../figure_code/'); % phi-three loading function is here

results_location = 'results/';
results_file = ['phi3_svm_' class_type '_singleFeature.mat'];

%% Load phi values

if strcmp(class_type, 'across')
    global_tpm = 1;
else % strcmp(class_type, 'within')
    global_tpm = 0;
end

[phi_threes, measure_strings{1}] = phi_load('phi_three', global_tpm, bin_location);

values_all = phi_threes{3}.phis(:, :, :, :, tau);
channel_sets = phi_threes{3}.channel_sets;
nSets = size(channel_sets, 1);

%% Plot

flies = [1 2 3 4 5 6 7 8 9 10 12 13];

figure;

xpos = zeros(length(flies), 1) + 1;

for set_counter = 1 : 100%size(channel_sets, 1)
    % Values
    scatter(xpos+0.1, squeeze(values_all(set_counter, 1, flies, 1)), 'r.'); hold on;
    scatter(xpos-0.1, squeeze(values_all(set_counter, 1, flies, 2)), 'b.');
    % Means
    scatter(xpos(1)+0.1, squeeze(mean(values_all(set_counter, 1, flies, 1), 3)), 'rs');
    scatter(xpos(1)-0.1, squeeze(mean(values_all(set_counter, 1, flies, 2), 3)), 'bs');
    
    xpos = xpos+1;
end

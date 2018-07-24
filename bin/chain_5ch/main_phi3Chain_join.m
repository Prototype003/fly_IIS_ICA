%% Description

%{

Joins results from phi-3 chaining

%}

%% Setup

flies = (1:13);
conditions = (1:2);
tau = 4;
trials = (1:8);
nChannels = 5; n_states = 2^nChannels;
global_tpm = 0;

bin_location = '../';

chain_type = 'phi'; % 'accuracy' for chaining using highest classification accuracy channels; 'phi' for using highest phi values
tau_type = 'tauBin'; % 'tau' if tau represents number of sample-steps; 'tauBin' if tau represents average over tau samples

if strcmp(chain_type, 'accuracy')
    source_directory = ['results_phiChain_split_accuracy/'];
    source_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels' num2str(nChannels) '_globalTPM' num2str(global_tpm) '_chainedPhi'];
    source_suffix = '.mat';
    
    output_directory = ['results_phiChain/'];
    if ~isdir(output_directory)
        mkdir(output_directory);
    end
    output_file = [source_prefix '_accuracy' source_suffix];
else % strcmp(chain_type, 'phi')
    source_directory = ['results_phiChain_split/'];
    source_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels' num2str(nChannels) '_globalTPM' num2str(global_tpm) '_chainedPhi'];
    source_suffix = '.mat';
    
    output_directory = ['results_phiChain/'];
    if ~isdir(output_directory)
        mkdir(output_directory);
    end
    output_file = [source_prefix '_' tau_type num2str(tau) source_suffix];
end
%% Join flies

% Create storage structures
channel_sets = zeros(nChannels, length(trials), length(flies));
phi_joined = zeros(length(trials), length(flies), length(conditions));
state_counters = zeros(n_states, length(trials), length(flies), length(conditions));
state_phis = zeros(n_states, length(trials), length(flies), length(conditions));
tpms = zeros(n_states, n_states, length(trials), length(flies), length(conditions));
mips = cell(n_states, length(trials), length(flies), length(conditions));

for fly = 1 : length(flies)
    disp(fly);
    for trial = 1 : length(trials)
        file_id = ['_f' num2str(fly, '%02d') tau_type num2str(tau) 't' num2str(trial)];
        
        % Load file
        load([source_directory source_prefix file_id source_suffix]);
        
        % Add to storage structures
        for condition = 1 : length(conditions)
            channel_sets(:, trial, fly) = phi.channel_set;
            phi_joined(trial, fly, condition) = phi.phis(condition);
            state_counters(:, trial, fly, condition) = phi.state_counters(:, condition);
            state_phis(:, trial, fly, condition) = phi.state_phis(:, condition);
            tpms(:, :, trial, fly, condition) = phi.tpms(:, :, condition);
            mips(:, trial, fly, condition) = phi.mips(:, condition);
        end
        
    end
end

phi_values = phi_joined;
phi_joined = struct();
phi_joined.channel_sets = channel_sets+1; % Fix python indexing
phi_joined.phis = phi_values;
phi_joined.state_counters = state_counters;
phi_joined.state_phis = state_phis;
phi_joined.tpms = tpms;
phi_joined.mips = mips;
phi_joined.nChannels = nChannels;
phi_joined.tau = tau;
phi_joined.global_tpm = 0;

%% Plot all trials

condition_colours = 'rb';
condition_offsets = [-0.05 0.05];

figure;
for condition = 1 : size(phi_joined.phis, 3)
    for trial = 1 : size(phi_joined.phis, 1)
        values = permute(phi_joined.phis(trial, :, condition), [2 1 3]);
        scatter(flies+condition_offsets(condition), values, 100, [condition_colours(condition) '.']); hold on;
    end
end

%% Plot trial averaged

condition_colours = 'rb';
values = permute(mean(phi_joined.phis, 1), [2 3 1]);

figure;
for condition = 1 : size(values, 2)
    scatter(flies, values(:, condition), condition_colours(condition)); hold on;
end


%% Compare to 2, 3, 4 channels

% addpath([bin_location 'figure_code/']);
% 
% [phis, measure_string] = phi_load('phi_three', 0, bin_location);
% 
% figure;
% 
% for nChannels_counter = 1 : 3
%     values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, :, :, 1), 2), 1), [3 4 1 2 5]); % Mean across all sets
%     %values = permute(max(mean(phis{nChannels_counter}.phis(:, :, :, :, 1), 2), [], 1), [3 4 1 2 5]); % Take the best sets (average across trials)
%     errorbar((1:2), mean(values, 1), std(values, [], 1) / sqrt(size(values, 1))); hold on;
% end
% 
% values = permute(mean(phi_joined.phis, 1), [2 3 1]);
% errorbar((1:2), mean(values, 1), std(values, [], 1) / sqrt(size(values, 1)));
% 
% axis_defaults(gca);
% legend('2ch', '3ch', '4ch', '5ch')
% xlim([0.6 2.4]);
% set(gca, 'XTick', [1 2], 'XTickLabel', {'W', 'A'});
% xlabel('condition');
% ylabel(measure_string);

%% Classification

accuracies = zeros(size(flies));

for fly = 1 : length(flies)
    
    class_data = permute(phi_joined.phis(:, fly, :), [1 3 2]);
    
    [accuracies(fly), ~, ~] = classify_lol(class_data);
    
end

phi_joined.accuracies = accuracies;

%% Save joined results

save([output_directory output_file], 'phi_joined');
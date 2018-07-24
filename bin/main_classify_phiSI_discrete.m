%% DESCRIPTION

%{

Classify awake/anest based on phi value

%}

%% Setup

flies = (1:13);
conditions = (1:2);
nChannels = (2:4);
tau = 1; % tau-index; 1=4ms; 2=8ms; 3=16ms

nFlies = length(flies);
nConditions = length(conditions);

data_detrended = 0;
data_zscored = 0;

deviation_type = 'nonpara';

global_covariance = 'nonGlobal'; % 'nonGlobal' 8 trials, or 'global' 1 big trial
if strcmp(global_covariance, 'global')
    global_tpm = 1;
    trials = 1;
    across_flies = 1;
else % strcmp(global_covariance, 'nonGlobal');
    global_tpm = 0;
    trials = (1:8);
    across_flies = 0;
end

addpath('figure_code/');

results_directory = 'workspace_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_phiSI_discrete_' global_covariance '_classification_across' num2str(across_flies) '.mat'];

%% Load
measure = 'phi_SI'; % 'phi_three' or 'phi_star'

[phis, measure_string] = phi_load(measure, global_tpm, '');

%% Classify

accuracies = cell(1, length(nChannels));

for nChannels_counter = 1 : length(nChannels)
    sets = (1 : size(phis{nChannels_counter}.channel_sets));
    nSets = length(sets);
    
    if across_flies == 1
        
        accuracy = zeros(nSets, 1);
        deviations = zeros(nSets, 1);
        
        for set_counter = 1 : nSets
            set = sets(set_counter);
            disp(set);
            
            % Classes are awake and anest, each column is a class
            class_data = [...
                permute(phis{nChannels_counter}.phis(set, trials, :, 1, tau), [3 1 2 4 5])...
                permute(phis{nChannels_counter}.phis(set, trials, :, 2, tau), [3 1 2 4 5])...
                ];
            
            % Classify
            [accuracy(set_counter), classifications, correctness] = classify_lol(class_data);
            
            % Error deviation
            deviations(set_counter) = error_distance(class_data, deviation_type);
        end
        
    else
        
        accuracy = zeros(nSets, nFlies);
        deviations = zeros(nSets, nFlies);
        
        for fly = 1 : nFlies
            for set_counter = 1 : nSets
                set = sets(set_counter);
                
                disp(['fly ' num2str(fly) ' set ' num2str(set)]);
                
                % Classes are awake and anest, each column is a class
                class_data = [...
                    permute(phis{nChannels_counter}.phis(set, trials, fly, 1, tau), [2 1 3 4 5])...
                    permute(phis{nChannels_counter}.phis(set, trials, fly, 2, tau), [2 1 3 4 5])...
                    ];
                
                % Classify
                [accuracy(set_counter, fly), classifications, correctness] = classify_lol(class_data);
                
                % Error deviation
                deviations(set_counter, fly) = error_distance(class_data, deviation_type);
            end
        end
        
    end
    
    accuracies{nChannels_counter}.nChannels = nChannels(nChannels_counter);
    accuracies{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
    accuracies{nChannels_counter}.accuracies = accuracy;
    accuracies{nChannels_counter}.deviations = deviations;
    
end

%% Plot

alpha = 0.25;
nChannels_counter = 3;

figure;
for fly = 1 : nFlies
    subplot(4, 4, fly);
    for trial = 1 : length(trials)
        scatter((1:size(phis{nChannels_counter}.phis, 1)), phis{nChannels_counter}.phis(:, trial, fly, 2, tau), 'b.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha); hold on;
        scatter((1:size(phis{nChannels_counter}.phis, 1)), phis{nChannels_counter}.phis(:, trial, fly, 1, tau), 'r.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    end
    title(['Fly ' num2str(fly)]);
end

% figure;
% plot(frequencies, accuracies(:, channel)); hold on;
% title('Classification (anest/awake) using power');
% xlabel('frequency (Hz)');
% ylabel('% correct');

%% Save

disp('saving');
save([results_directory results_filename], 'accuracies');
disp('saved');
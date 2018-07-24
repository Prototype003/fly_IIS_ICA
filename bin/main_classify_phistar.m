%% DESCRIPTION

%{

Classify awake/anest based on phi value

%}

%% Setup

flies = (1:13);
conditions = (1:2);
nChannels = (2:4);
tau = 1;
deviation_type = 'nonpara'; % 'para' for sum of error deviations or 'nonpara' for rank correlation from ideal order

nFlies = length(flies);
nConditions = length(conditions);

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

global_covariance = 'global'; % 'nonGlobal' Eight trials, each with own covariance matrices, or 'global' 1 covariance matrix built across all 8 trials

if strcmp(global_covariance, 'nonGlobal')
    data_directory = 'results/';
    data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
        '_detrend' num2str(data_detrended)...
        '_zscore' num2str(data_zscored)...
        '_nChannels' data_nChannels...
        %'_shareFiltered'
        ];
    data_filename = [data_filename '_phistar.mat'];
    trials = (1:8); % 8 trials
    across_flies = 0;
else % strcmp(global_covariance, 'global')
    data_directory = 'results/';
    data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_medianSplit0_phistar_global';
    data_filename = [data_filename '.mat'];
    trials = 1; % Global calculation means only 1 effective trial (identical covariance across 8 trials)
    across_flies = 1;
end
results_directory = 'workspace_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_phistar_' global_covariance '_classification_across' num2str(across_flies) '.mat'];

%% Load

disp('Loading');
load([data_directory data_filename]);
disp('Loaded');

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
                permute(phis{nChannels_counter}.phi_stars(set, trials, :, 1, tau), [3 1 2 4 5])...
                permute(phis{nChannels_counter}.phi_stars(set, trials, :, 2, tau), [3 1 2 4 5])...
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
                    permute(phis{nChannels_counter}.phi_stars(set, trials, fly, 1, tau), [2 1 3 4 5])...
                    permute(phis{nChannels_counter}.phi_stars(set, trials, fly, 2, tau), [2 1 3 4 5])...
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
        scatter((1:size(phis{nChannels_counter}.phi_stars, 1)), phis{nChannels_counter}.phi_stars(:, trial, fly, 2, tau), 'b.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha); hold on;
        scatter((1:size(phis{nChannels_counter}.phi_stars, 1)), phis{nChannels_counter}.phi_stars(:, trial, fly, 1, tau), 'r.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
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
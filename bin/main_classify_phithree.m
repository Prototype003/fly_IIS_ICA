%% DESCRIPTION

%{

Classify awake/anest based on phi value

%}

%% Setup

flies = (1:13);
conditions = (1:2);
nChannels = (2:4);
tau = 1;
deviation_type = 'nonpara'; % 'para' for sum of error deviations or 'nonpara' for rank correlation

nFlies = length(flies);
nConditions = length(conditions);

data_nChannels = {'2t2', '3t3', '4t4'}; %, '3t3', '4t4'};
data_detrended = 0;
data_zscored = 0;

across_flies = 0;
global_tpm = 'global';

results_directory = 'workspace_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_phithree_' global_tpm '_classification_across' num2str(across_flies) '.mat'];

%% Classify

accuracies = cell(1, length(nChannels));
for nChannels_counter = 1 : length(data_nChannels)
    
    nChannels_infix = data_nChannels{nChannels_counter};
    data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
        '_detrend' num2str(data_detrended)...
        '_zscore' num2str(data_zscored)...
        '_nChannels' nChannels_infix...
        %'_shareFiltered'
        ];
    
    % Load file
    disp('loading')
    if strcmp(global_tpm, 'nonGlobal')
        if strcmp(nChannels_infix, '2t2')
            data_directory = 'results/';
            data_filename = [data_filename '_phithree_nonGlobal.mat'];
            load([data_directory data_filename]);
            phis = phis{1};
        elseif strcmp(nChannels_infix, '3t3')
            data_directory = 'results/';
            data_filename = [data_filename '_phithree_nonGlobal_tau4.mat'];
            load([data_directory data_filename]);
            phis = phis{1};
        else % '4t4'
            data_directory = 'results_split/';
            data_filename = [data_filename '_phithree_nonGlobal_tau4.mat'];
            phis = load([data_directory data_filename]);
            phis.channel_sets = nchoosek((1:15), 4);
        end
    else
        data_directory = 'results/';
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree.mat'];
        load([data_directory data_filename]);
        phis = phis{nChannels_counter};
    end
    disp('loaded');
    
    sets = (1 : size(phis.channel_sets));
    nSets = length(sets);
    
    % Classify
    
    if across_flies == 1
        
        accuracy = zeros(nSets, 1);
        deviations = zeros(nSets, 1);
        
        % Average across trials to get a single phi value
        phis.phi_threes = mean(phis.phi_threes, 2);
        
        for set_counter = 1 : nSets
            set = sets(set_counter);
            disp(['set ' num2str(set)]);
            
            % Classes are awake and anest, each column is a class
            class_data = [...
                permute(phis.phi_threes(set, 1, :, 1, tau), [3 1 2 4 5])...
                permute(phis.phi_threes(set, 1, :, 2, tau), [3 1 2 4 5])...
                ];
            
            % Classify
            [accuracy(set_counter), leave_out, classifications, correctness] = classify_lol(class_data);
            
            % Error deviation
            deviations(set_counter) = error_distance(class_data, deviation_type);
        end
        
    else % within flies
        
        accuracy = zeros(nSets, nFlies);
        deviations = zeros(nSets, nFlies);
        
        for fly = 1 : nFlies
            for set_counter = 1 : nSets
                set = sets(set_counter);
                
                disp(['fly ' num2str(fly) ' set ' num2str(set)]);
                
                % Classes are awake and anest, each column is a class
                class_data = [...
                    permute(phis.phi_threes(set, :, fly, 1, tau), [2 1 3 4 5])...
                    permute(phis.phi_threes(set, :, fly, 2, tau), [2 1 3 4 5])...
                    ];
                
                % Classify
                [accuracy(set_counter, fly), leave_out, classifications, correctness] = classify_lol(class_data);
                
                % Error deviation
                deviations(set_counter, fly) = error_distance(class_data, deviation_type);
            end
        end
    end
    
    accuracies{nChannels_counter}.nChannels = nChannels(nChannels_counter);
    accuracies{nChannels_counter}.channel_sets = phis.channel_sets;
    accuracies{nChannels_counter}.accuracies = accuracy;
    accuracies{nChannels_counter}.deviations = deviations;
    
end

%% Plot
fly = 8;
alpha = 0.25;

figure;
for fly = 1 : nFlies
    subplot(4, 4, fly);
    for trial = 1 : size(phis.phi_threes, 2)
        scatter((1:size(phis.phi_threes, 1)), phis.phi_threes(:, trial, fly, 2, tau), 'b.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha); hold on;
        scatter((1:size(phis.phi_threes, 1)), phis.phi_threes(:, trial, fly, 1, tau), 'r.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    end
    title(['Fly ' num2str(fly)]);
end

%% Average across flies
figure;
for trial = 1 : size(phis.phi_threes, 2)
    scatter((1:size(phis.phi_threes, 1)), mean(phis.phi_threes(:, trial, :, 2, tau), 3), 'b.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha); hold on;
    scatter((1:size(phis.phi_threes, 1)), mean(phis.phi_threes(:, trial, :, 1, tau), 3), 'r.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
end

figure;
scatter((1:size(phis.phi_threes, 1)), mean(mean(phis.phi_threes(:, :, :, 1, tau), 2)-mean(phis.phi_threes(:, :, :, 2, tau), 2), 3), 'b.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha); hold on;
%scatter((1:size(phis.phi_threes, 1)), mean(mean(phis.phi_threes(:, :, fly, 1, tau), 2), 3), 'r.', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);


%% Save

disp('Saving');
save([results_directory results_filename], 'accuracies');
disp('Saved');

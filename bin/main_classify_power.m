%% DESCRIPTION

%{

Calculate power for single channels and classify awake/anest based on power

%}

%% Setup

flies = (1:13);
conditions = (1:2);
channels = (1:15);

chronux_params = struct();
chronux_params.tapers = [2 3]; % 3 tapers, 3 = 2*2 - 1
chronux_params.Fs = 1000; % 1000 Hz downsampled sampling rate of the input data
chronux_params.fpass = [0 50]; % This is the frequency range we're interested in
chronux_params.pad = 1; % I think Dror used 1, apparently higher padding gives higher frequency resolution (more frequency bins) but doesn't affect calculation
% .err - default no errorbars
% .trialave - default no averaging

nFlies = length(flies);
nConditions = length(conditions);
nChannels = length(channels);

across_flies = 1;

data_directory = 'workspace_results/';
data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat';

results_directory = 'workspace_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across' num2str(across_flies) '.mat'];

%% Load

disp('loading');
load([data_directory data_filename]);
disp('loaded');

%% Calculate power
% Chronux continuous data - first dimension is time, second dimension is trials/channels

disp('Calculating power');

if across_flies == 0
    
    % Results matrix: frequencies x trials x channels x condition x fly
    powers = zeros(410, size(fly_data, 3), nChannels, nConditions, nFlies);
    
    for condition = 1 : nConditions
        for fly_counter = 1 : nFlies
            fly = flies(fly_counter);
            
            for channel_counter = 1 : nChannels
                channel = channels(channel_counter);
                
                % Get trials for the channel (input dimensions to chronux function should be time x trials
                trials = permute(fly_data(:, channel, :, fly, condition), [1 3 2 4 5]);
                
                % Compute power spectrums
                [spectrums, frequencies] = mtspectrumc(trials, chronux_params);
                
                % Store
                powers(:, :, channel_counter, condition, fly_counter) = spectrums;
            end
            
        end
    end
    
else % across == 1
    
    % Results matrix: frequencies x trials x channels x condition x fly
    powers = zeros(3277, 1, nChannels, nConditions, nFlies);
    
    for condition = 1 : nConditions
        for fly_counter = 1 : nFlies
            fly = flies(fly_counter);
            
            for channel_counter = 1 : nChannels
                channel = channels(channel_counter);
                
                % Get trials for the channel (input dimensions to chronux function should be time x trials
                trials = permute(fly_data(:, channel, :, fly, condition), [1 3 2 4 5]);
                
                % Concatenate trials into one big trial
                trial = zeros(size(trials, 1) * size(trials, 2), size(trials, 3));
                row_counter = (1:size(trials, 1));
                for trial_counter = 1 : size(trials, 2)
                    trial(row_counter, :) = trials(:, trial_counter, :);
                    row_counter = row_counter + size(trials, 1);
                end
                
                % Compute power spectrums
                [spectrums, frequencies] = mtspectrumc(trial, chronux_params);
                
                % Store
                powers(:, :, channel_counter, condition, fly_counter) = spectrums;
            end
            
        end
    end
    
end

powers = log(powers);

disp('log powers calculated');

%% Classify

disp('classifying');

if across_flies == 1 % across flies
    
    powers = mean(powers, 2);
    
    accuracies = zeros(length(frequencies), nChannels);
    classifications = cell(1, nChannels);
    
    for channel = 1 : nChannels
        disp(['channel ' num2str(channel)]);
        for frequency = 1 : length(frequencies)
            
            % Classes are awake and anest, each column is a class
            class_data = [...
                permute(powers(frequency, 1, channel, 1, :), [5 1 2 3 4])...
                permute(powers(frequency, 1, channel, 2, :), [5 1 2 3 4])...
                ];
            [accuracies(frequency, channel), classifications{channel}.leave_out, classifications{channel}.classifications, classifications{channel}.correctness] = classify_lol(class_data);
        end
    end
    
else % within flies
    
    accuracies = zeros(length(frequencies), nChannels, nFlies);
    classifications = cell(1, nFlies);
    
    for fly = 1 : nFlies
        for channel = 1 : nChannels
            disp(['fly ' num2str(fly) ' channel ' num2str(channel)]);
            for frequency = 1 : length(frequencies)
                
                % Classes are awake and anest, each column is a class
                class_data = [...
                    permute(powers(frequency, :, channel, 1, fly), [2 1 3 4 5])...
                    permute(powers(frequency, :, channel, 2, fly), [2 1 3 4 5])...
                    ];
                
                [accuracies(frequency, channel, fly), classifications{fly}.leave_out, classifications{fly}.classifications, classifications{fly}.correctness] = classify_lol(class_data);
            end
        end
    end
    
end

disp('classified');

%% Save

save([results_directory results_filename], 'chronux_params', 'powers', 'frequencies', 'accuracies', 'classifications');

%% Plot average

figure;
errorbar(...
    frequencies,...
    mean(mean(accuracies, 2), 3),...
    0.5*std(mean(accuracies, 2), [], 3) / sqrt(nFlies));

%% Plot for 1 fly

fly = 3;
figure;
for channel = 1 : nChannels
    plot(frequencies, accuracies(:, channel, fly)); hold on;
end
title('Classification (anest/awake) using power');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass 0 100]);

% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);

%% Plot for 1 channel (for across fly classification)

channel = 1;

figure;
for channel = 1 : nChannels
    subplot(4, 4, channel);
    for fly = 1 : 13
        scatter((1:410), squeeze(powers(:, 1, channel, 1, fly)), 'r.'); hold on;
        scatter((1:410), squeeze(powers(:, 1, channel, 2, fly)), 'b.');
    end
end

%%
figure;
for channel = 1 : nChannels
    plot((1:410), accuracies(:, channel)); hold on;
end
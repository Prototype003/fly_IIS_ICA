%% DESCRIPTION

%{

Calculate power for single channels and classify awake/anest based on power

%}

%% Setup

flies = (1:13);
conditions = (1:2);
nChannels = 15;
network_size = 2;

chronux_params = struct();
chronux_params.tapers = [3 5]; % 5 tapers, 3 = 2*3 - 1
chronux_params.Fs = 1000; % 1000 Hz downsampled sampling rate of the input data
chronux_params.fpass = [0 50]; % This is the frequency range we're interested in
chronux_params.pad = 1; % I think Dror used 1, apparently higher padding gives higher frequency resolution (more frequency bins) but doesn't affect calculation
% .err - default no errorbars
% .trialave - default no averaging

nFlies = length(flies);
nConditions = length(conditions);

across_flies = 1;

data_directory = 'workspace_results/';
data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat';

results_directory = 'workspace_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification_across' num2str(across_flies) '.mat'];

%% Load

load([data_directory data_filename]);

%% Calculate coherence
disp('Computing coherencies');

% Chronux continuous data - first dimension is time, second dimension is trials/channels
networks = nchoosek((1:nChannels), network_size);
nNetworks = size(networks, 1);

% Results matrix: frequencies x trials x channel sets x condition x fly
coherencies = zeros(410, size(fly_data, 3), nNetworks, nConditions, nFlies);

for condition = 1 : nConditions
    for fly_counter = 1 : nFlies
        fly = flies(fly_counter);
        
        for network_counter = 1 : nNetworks
            network = networks(network_counter, :);
            disp(network);
            
            % Get trials for the channels in the network (input dimensions to chronux function should be time x trials
            trials = permute(fly_data(:, network, :, fly, condition), [1 3 2 4 5]);
            
            % Compute power spectrums
            [spectrums, ~, ~, ~, ~, frequencies] = coherencyc(trials(:, :, 1), trials(:, :, 2), chronux_params);
            
            % Store
            coherencies(:, :, network_counter, condition, fly_counter) = spectrums;
        end
        
    end
end

disp('Coherencies computed');

%% Classify
disp('Classifying');

if across_flies == 1
    
    corherencies = mean(coherencies, 2);
    
    accuracies = zeros(length(frequencies), nNetworks);
    classifications = cell(1, nNetworks);
    
    for network = 1 : nNetworks
        disp(['Pair ' num2str(network)]);
        for frequency = 1 : length(frequencies)
            
            % Classes are awake and anest, each column is a class
            class_data = [...
                permute(coherencies(frequency, 1, network, 1, :), [5 1 2 3 4])...
                permute(coherencies(frequency, 1, network, 2, :), [5 1 2 3 4])...
                ];
            
            [accuracies(frequency, network), classifications{network}.leave_out, classifications{network}.classifications, classifications{network}.correctness] = classify_lol(class_data);
            
        end
    end
    
else % within flies
    
    accuracies = zeros(length(frequencies), nNetworks, nFlies);
    classifications = cell(1, nFlies);
    
    for fly = 1 : nFlies
        for network = 1 : nNetworks
            disp(['Fly ' num2str(fly) ' network ' num2str(network)]);
            for frequency = 1 : length(frequencies)
                
                % Classes are awake and anest, each column is a class
                class_data = [...
                    permute(coherencies(frequency, :, network, 1, fly), [2 1 3 4 5])...
                    permute(coherencies(frequency, :, network, 2, fly), [2 1 3 4 5])...
                    ];
                
                [accuracies(frequency, network, fly), classifications{fly}.leave_out, classifications{fly}.classifications, classifications{fly}.correctness] = classify_lol(class_data);
            end
        end
    end
    
end

disp('Classified');

%% Save

save([results_directory results_filename], 'chronux_params', 'coherencies', 'networks', 'frequencies', 'accuracies', 'classifications');

%% Plot

figure;
for network = 1 : nNetworks
    plot(frequencies, accuracies(:, network)); hold on;
end
title('Classification (anest/awake) using coherence');
xlabel('frequency (Hz)');
ylabel('% correct');

%% Description

%{

Classifies wake/anest using power, coherence, with parameters to match phi

Power: use n features (power of each channel) for n channels

Coherence: use n features (pairwise coherence between each pair of channels)
    for n channels, pairwise features correspond to 2nd order concepts

Phi - classification using all concepts is already conducted in script
    ../phi_3/main_composition_classify.m


%}

%% SETUP

nChannels = (2:4);

class_type = 'across'; % 'across' or 'within'

freq_range = (1:42); % (1:42) = 0-5Hz; Vector of which frequency bins to use, we want to use the lowest frequency bin
%freq_range = (2:42); % (2:42) = 0.1-5Hz;
%freq_range = (1:83); % (1:83) = 0-10Hz;
%freq_range = (1:410); % Everything

freq_range = (83:165);

% bin directory location
bin_dir = '../';

results_location = 'results/';
results_file = ['power_svm_channelMatched_' class_type '.mat'];

%% Load power

powers = load_power(bin_dir);

% Mean across frequency bins
powers = mean(powers(freq_range, :, :, :, :), 1);
powers = permute(powers, [2 3 4 5 1]); % trials x channels x conditions x flies

%% Load without average across frequency range
powers = load_power(bin_dir);

%% Setup features + results storage structure

power_accuracies = cell(length(nChannels), 1);

for nChannels_counter = 1 : length(nChannels)
    power_accuracies{nChannels_counter}.nChannels = nChannels(nChannels_counter);
    power_accuracies{nChannels_counter}.channel_sets = nchoosek((1:size(powers, 2)), nChannels(nChannels_counter));
end

%% Classify (within flies)

% Uses SVM

cost = 1;

for nChannels_counter = 1 : length(nChannels)
    
    % Results storage
    accuracies = zeros(...
        length(freq_range),...
        size(powers, 4),...
        size(power_accuracies{nChannels_counter}.channel_sets, 1)...
        );
    
    for freq = 1 : length(freq_range)
        disp([num2str(nChannels_counter) ' ' num2str(freq)]); tic;
        for fly = 1 : size(powers, 5)
            
            for network = 1 : size(power_accuracies{nChannels_counter}.channel_sets,1)
                
                channel_set = power_accuracies{nChannels_counter}.channel_sets(network, :);
                
                % Get features
                features = permute(powers(freq, :, channel_set, :, fly), [2 3 4 1 5]); % All features
                %features = permute(mean(powers(freq, :, channel_set, :, fly), 3), [2 4 3 1]); % Mean across channels
                
                % Classify
                results = svm_lol_liblinear_manual(features, cost);
                accuracies(freq, fly, network) = results.accuracy;
                
            end
            
        end
        toc
    end
    
    power_accuracies{nChannels_counter}.accuracies = accuracies;
    
end

%% Plot;


figure;

%% Classify (within flies)

% Uses nearest-mean classification, on averaged power across channels

addpath(bin_dir);

for nChannels_counter = 1 : length(nChannels)
    
    % Results storage
    accuracies = zeros(...
        length(freq_range),...
        size(powers, 4),...
        size(power_accuracies{nChannels_counter}.channel_sets, 1)...
        );
    
    for freq_counter = 1 : length(freq_range)
        freq = freq_range(freq_counter);
        disp([num2str(nChannels_counter) ' ' num2str(freq)]); tic;
        for fly = 1 : size(powers, 5)
            
            for network = 1 : size(power_accuracies{nChannels_counter}.channel_sets,1)
                
                channel_set = power_accuracies{nChannels_counter}.channel_sets(network, :);
                
                % Get features
                features = permute(mean(powers(freq, :, channel_set, :, fly), 3), [2 4 3 1]);
                
                % Classify
                results = classify_lol(features);
                accuracies(freq_counter, fly, network) = results;
                
            end
            
        end
        toc
    end
    
    power_accuracies{nChannels_counter}.accuracies = accuracies;
    
end

%save('power_nChannelsMeaned_classify_across0.mat', 'power_accuracies', 'nChannels');

%% Classify (across flies)

cost = 1;

for nChannels_counter = 1 : length(nChannels)
    
    % Results storage
    accuracies = zeros(size(power_accuracies{nChannels_counter}.channel_sets, 1), 1);
    
    disp(num2str(nChannels_counter)); tic;
    
    for network = 1 : size(power_accuracies{nChannels_counter}.channel_sets,1)
        
        channel_set = power_accuracies{nChannels_counter}.channel_sets(network, :);
        
        % Get features
        features = permute(mean(powers(:, channel_set, :, :), 1), [4 2 3 1]);
        features = permute(mean(mean(powers(:, channel_set, :, :), 1), 2), [4 2 3 1]);
        
        % Classify
        results = svm_lol_liblinear_manual(features, cost);
        accuracies(network) = results.accuracy;
        
    end
    toc
    power_accuracies{nChannels_counter}.accuracies = accuracies;
    
end

%% Classify (across flies)

% Uses nearest-mean classification, on averaged power across channels

for nChannels_counter = 1 : length(nChannels)
    
    % Results storage
    accuracies = zeros(length(freq_range), size(power_accuracies{nChannels_counter}.channel_sets, 1));
    
    for freq_counter = 1 : length(freq_range)
        freq = freq_range(freq_counter);
        disp([num2str(nChannels_counter) ' ' num2str(freq)]); tic;
        
        for network = 1 : size(power_accuracies{nChannels_counter}.channel_sets,1)
            
            channel_set = power_accuracies{nChannels_counter}.channel_sets(network, :);
            
            % Get features
            features = permute(mean(mean(powers(freq, :, channel_set, :, :), 3), 2), [5 4 3 2 1]);
            
            % Classify using nearest-mean
            %results = classify_lol(features);
            %accuracies(freq_counter, network) = results;
            
            % Classify using svm
            results = svm_lol_liblinear_manual(permute(features, [1 3 2]), cost);
            accuracies(freq_counter, network) = results.accuracy;
            
        end
        toc
    end
    power_accuracies{nChannels_counter}.accuracies = accuracies;
    
end

save('power_nChannelsMeaned_classify_across1.mat', 'power_accuracies', 'nChannels', 'freq_range');

%% Save accuracies

results_dir = 'results/';
results_file = ['power_classify_across' num2str(across) '.mat'];
save([results_dir results_file], 'power_accuracies', 'nChannels');
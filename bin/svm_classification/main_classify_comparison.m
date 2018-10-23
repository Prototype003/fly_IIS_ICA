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
freq_range = (2:42); % (2:42) = 0.1-5Hz;
freq_range = (1:83); % (1:83) = 0-10Hz;
freq_range = (1:410); % Everything

% bin directory location
bin_dir = '../';

results_location = 'results/';
results_file = ['power_svm_channelMatched_' class_type '.mat'];

%% Load power

powers = load_power(bin_dir);

% Mean across frequency bins
powers = mean(powers(freq_range, :, :, :, :), 1);
powers = permute(powers, [2 3 4 5 1]); % trials x channels x conditions x flies

%% Setup features + results storage structure

power_accuracies = cell(length(nChannels), 1);

for nChannels_counter = 1 : length(nChannels)
    power_accuracies{nChannels_counter}.nChannels = nChannels(nChannels_counter);
    power_accuracies{nChannels_counter}.channel_sets = nchoosek((1:size(powers, 2)), nChannels(nChannels_counter));
end

%% Classify (within flies)

cost = 1;

for nChannels_counter = 1 : length(nChannels)
    
    % Results storage
    accuracies = zeros(...
        size(powers, 4),...
        size(power_accuracies{nChannels_counter}.channel_sets, 1)...
        );
    
    for fly = 1 : size(powers, 4)
        disp([num2str(nChannels_counter) ' ' num2str(fly)]); tic;
        for network = 1 : size(power_accuracies{nChannels_counter}.channel_sets,1)
            
            channel_set = power_accuracies{nChannels_counter}.channel_sets(network, :);
            
            % Get features
            features = powers(:, channel_set, :, fly);
            
            % Classify
            results = svm_lol_liblinear_manual(features, cost);
            accuracies(fly, network) = results.accuracy;
            
        end
        toc
    end
    
    power_accuracies{nChannels_counter}.accuracies = accuracies;
    
end

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
        
        % Classify
        results = svm_lol_liblinear_manual(features, cost);
        accuracies(network) = results.accuracy;
        
    end
    toc
    power_accuracies{nChannels_counter}.accuracies = accuracies;
    
end

%% Save accuracies

results_dir = 'results/';
results_file = ['composition_' constellation_type '_classify.mat'];
save([results_dir results_file], 'accuracies', 'tau', 'nChannels');
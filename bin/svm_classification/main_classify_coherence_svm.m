%% DESCRIPTION

%{

Classifies awake/anest using 105 features - cohereence of 15 channels (at lowest frequency)

Across flies classification
Within flies classification

%}

%% SETUP

class_type = 'within'; % 'across' or 'within'

freq_range = (1:42); % (1:42) = 0-5Hz; Vector of which frequency bins to use, we want to use the lowest frequency bin
freq_range = (43:165); % 10-20Hz

% bin directory location
bin_dir = '../';

results_location = 'results/';
results_file = ['coherence_svm_' class_type '.mat'];

%% Load power

coherencies = load_coherence(bin_dir);

% Mean across frequency bins
coherencies = mean(coherencies(freq_range, :, :, :, :), 1);
coherencies = permute(coherencies, [2 3 4 5 1]); % trials x channels x conditions x flies

%% Classify using SVM

if strcmp(class_type, 'across')
    values = permute(mean(coherencies, 1), [4 2 3 1]); % flies x channels x conditions
    
    classifications = svm_lol_libsvm(values);
    
elseif strcmp(class_type, 'within')
    
    classifications = cell(size(coherencies, 4), 1);
    for fly = 1 : size(coherencies, 4)
        disp(num2str(fly));
        values = coherencies(:, :, :, fly); % trials x channels x conditions
        
        classifications{fly} = svm_lol_libsvm(values);
    end
    
    accuracy_mean = 0;
    accuracy_flies = zeros(size(coherencies, 4), 1);
    for fly = 1 : size(coherencies, 4)
        accuracy_mean = accuracy_mean + classifications{fly}.accuracy;
        accuracy_flies(fly) = classifications{fly}.accuracy;
    end
    accuracy_mean = accuracy_mean / size(coherencies, 4)
    
end

%% Save

save([results_location results_file], 'classifications');
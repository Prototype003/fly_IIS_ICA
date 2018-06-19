%% DESCRIPTION

%{

Classifies awake/anest using 15 features - power of 15 channels (at lowest frequency)

Across flies classification
Within flies classification

%}

%% SETUP

class_type = 'within'; % 'across' or 'within'

freq_range = (1); % (1:42) = 0-5Hz; Vector of which frequency bins to use, we want to use the lowest frequency bin

% bin directory location
bin_dir = '../';

results_location = 'results/';
results_file = ['power_svm_' class_type '.mat'];

%% Load power

powers = load_power(bin_dir);

% Mean across frequency bins
powers = mean(powers(freq_range, :, :, :, :), 1);
powers = permute(powers, [2 3 4 5 1]); % trials x channels x conditions x flies

%% Classify using SVM

if strcmp(class_type, 'across')
    values = permute(mean(powers, 1), [4 2 3 1]); % flies x channels x conditions
    
    classifications = svm_lol_libsvm(values);
    
elseif strcmp(class_type, 'within')
    
    classifications = cell(size(powers, 4), 1);
    for fly = 1 : size(powers, 4)
        disp(num2str(fly));
        values = powers(:, :, :, fly); % trials x channels x conditions
        
        classifications{fly} = svm_lol_libsvm(values);
    end
    
end

%% Save

save([results_location results_file], 'classifications');
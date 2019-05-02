%% DESCRIPTION

%{

Classifies awake/anest using 105 features - correlation of 15 channels

Across flies classification
Within flies classification

%}

%% SETUP

class_type = 'across'; % 'across' or 'within'

% bin directory location
bin_dir = '../';

results_location = 'results/';
results_file = ['correlation_svm_' class_type '.mat'];

%% Compute correlations
% Matrix should be trials x pairs x conditions x flies

load('../workspace_results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat');

nTrials = size(fly_data, 3);
nFlies = size(fly_data, 4);
nConditions = size(fly_data, 5);
pairs = nchoosek((1:size(fly_data, 2)), 2);

% Can add pre-processing (e.g. binarisation) here
% Median split
fly_medians = median(fly_data, 1);
fly_medians = repmat(fly_medians, [size(fly_data, 1) 1 1 1 1]);
fly_data = fly_data > fly_medians;

% Correlation between each pair of channels
upper_triangle = logical(tril(ones(size(fly_data,2), size(fly_data,2)), -1)); % tril to match order of nchoosek (because indexing is row-wise first
correlations = zeros(nTrials, size(pairs, 1), nConditions, nFlies);
for fly = 1 : size(fly_data, 4)
    for condition = 1 : size(fly_data, 5)
        for trial = 1 : size(fly_data, 3)
            r = corr(fly_data(:, :, trial, fly, condition));
            correlations(trial, :, condition, fly) = r(upper_triangle);
        end
    end
end

%% Classify using SVM

if strcmp(class_type, 'across')
    values = permute(mean(correlations, 1), [4 2 3 1]); % flies x channels x conditions
    
    classifications = svm_lol_libsvm(values);
    
elseif strcmp(class_type, 'within')
    
    classifications = cell(size(correlations, 4), 1);
    for fly = 1 : size(correlations, 4)
        disp(num2str(fly));
        values = correlations(:, :, :, fly); % trials x channels x conditions
        
        classifications{fly} = svm_lol_libsvm(values);
    end
    
    accuracy_mean = 0;
    accuracy_flies = zeros(size(correlations, 4), 1);
    for fly = 1 : size(correlations, 4)
        accuracy_mean = accuracy_mean + classifications{fly}.accuracy;
        accuracy_flies(fly) = classifications{fly}.accuracy;
    end
    accuracy_mean = accuracy_mean / size(correlations, 4)
    
end

%% Save

save([results_location results_file], 'classifications');
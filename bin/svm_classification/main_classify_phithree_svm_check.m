%% DESCRIPTION

%{

Classifies awake/anest using 1 feature, similar to nearest-mean classification

Across flies classification
Within flies classification

%}

%% SETUP

class_type = 'within'; % 'across' or 'within'

nFeatures = 1;

tau = 1;

% bin directory location
bin_location = '../';

addpath('../figure_code/'); % phi-three loading function is here

results_location = 'results/';
results_file = ['phi3_svm_' class_type '_singleFeature.mat'];

%% Load phi values

if strcmp(class_type, 'across')
    global_tpm = 1;
else % strcmp(class_type, 'within')
    global_tpm = 0;
end

[phi_threes, measure_strings{1}] = phi_load('phi_three', global_tpm, bin_location);

values_all = phi_threes{3}.phis(:, :, :, :, tau);
channel_sets = phi_threes{3}.channel_sets;
nSets = size(channel_sets, 1);

%% Classify using SVM

if strcmp(class_type, 'across')
    
    classifications = cell(size(values_all, 1), 1);
    accuracy = zeros(nSets, 1);
    
    for set_counter = 1 : nSets
        disp(['set ' num2str(set_counter)]);
        
        values = permute(mean(values_all(set_counter, :, :, :), 2), [3 1 4 2]); % flies x set x conditions
        
        % Classify
        classifications{set_counter} = svm_lol_libsvm(values);
        
        % Summarise across validations and store
        classifications{set_counter}.accuracy_per_class = classifications{set_counter}.correct_total / size(classifications{set_counter}.leave_outs, 1);
        classifications{set_counter}.accuracy = sum(classifications{set_counter}.correct_total) / numel(classifications{set_counter}.leave_outs);
    end
    
elseif strcmp(class_type, 'within')
    
    classifications = cell(size(values_all, 3), 1);
    for fly = 1 : size(values_all, 3)
        disp(num2str(fly));
        
        % First classification round
        selected_sets = set_sample(channel_sets, nFeatures); % Select random set of channel sets
        values = permute(values_all(selected_sets, :, fly, :), [2 1 4 3]); % trials x sets x conditions
        classifications{fly} = svm_lol(values);
        
        % Store selected sets
        classifications{fly}.sets = zeros(size(classifications{fly}.leave_outs, 1), length(selected_sets));
        classifications{fly}.sets(1, :) = selected_sets;
        classifications{fly}.channel_sets = channel_sets;
        
        % Reset total counters
        classifications{fly}.correct_total = classifications{fly}.correctness(1, :);
        
        % Subsequent classification rounds
        for validation = 2 : size(classifications{fly}.leave_outs, 1)
            selected_sets = set_sample(channel_sets, nFeatures);
            classifications{fly}.sets(validation, :) = selected_sets;
            values = permute(values_all(selected_sets, :, fly, :), [2 1 4 3]); % trials x sets x conditions
            tmp = svm_lol(values, validation);
            
            % Store values used and results
            classifications{fly}.class_data(validation, :) = tmp.class_data;
            classifications{fly}.svms{validation} = tmp.svms{validation};
            classifications{fly}.predictions(validation, :) = tmp.predictions(validation, :);
            classifications{fly}.confidences(validation, :, :) = tmp.confidences(validation, :, :);
            classifications{fly}.correctness(validation, :) = tmp.correctness(validation, :);
            classifications{fly}.correct_total = classifications{fly}.correct_total + tmp.correctness(validation, :);
        end
        
        % Summarise across validations
        classifications{fly}.accuracy_per_class = classifications{fly}.correct_total / size(classifications{fly}.leave_outs, 1);
        classifications{fly}.accuracy = sum(classifications{fly}.correct_total) / numel(classifications{fly}.leave_outs);
    end
    
end

%% Save

save([results_location results_file], 'classifications');

%%

figure;

for fly = 1 : size(values_all, 3)
    subplot(size(values_all, 3), 1, fly);
    
    bar(squeeze(values_all(:, :, fly, 1) - values_all(:, :, fly, 2)));
end


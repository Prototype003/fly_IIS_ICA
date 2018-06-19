%% DESCRIPTION

%{

Classifies awake/anest using 15 features - 15 phis for 15 sets of 4 channels (spanning 15 channels)

Across flies classification
Within flies classification

%}

%% SETUP

class_type = 'within'; % 'across' or 'within'

nFeatures = 15;

tau = 1;

% bin directory location
bin_location = '../';

addpath('../figure_code/'); % phi-three loading function is here

results_location = 'results/';
results_file = ['phi3_svm_' class_type '.mat'];

%% Load phi values

if strcmp(class_type, 'across')
    global_tpm = 1;
else % strcmp(class_type, 'within')
    global_tpm = 0;
end

[phi_threes, measure_strings{1}] = phi_load('phi_three', global_tpm, bin_location);

values_all = phi_threes{3}.phis(:, :, :, :, tau);
channel_sets = phi_threes{3}.channel_sets;

%% Classify using SVM

if strcmp(class_type, 'across')
    
    % First classification round
    selected_sets = set_sample(channel_sets, nFeatures); % Select random set of channel sets
    values = permute(mean(values_all(selected_sets, :, :, :), 2), [3 1 4 2]); % flies x sets x conditions
    classifications = svm_lol(values);
    
    % Store selected sets
    classifications.sets = zeros(size(classifications.leave_outs, 1), length(selected_sets));
    classifications.sets(1, :) = selected_sets;
    classifications.channel_sets = channel_sets;
    
    % Reset total counters
    classifications.correct_total = classifications.correctness(1, :);
    
    % Subsequent classification rounds
    for validation = 2 : size(classifications.leave_outs, 1)
        selected_sets = set_sample(channel_sets, nFeatures);
        classifications.sets(validation, :) = selected_sets;
        values = permute(mean(values_all(selected_sets, :, :, :), 2), [3 1 4 2]); % flies x sets x conditions
        tmp = svm_lol_libsvm(values, validation);
        
        % Store values used and results
        classifications.class_data(validation, :) = tmp.class_data;
        classifications.svms{validation} = tmp.svms{validation};
        classifications.predictions(validation, :) = tmp.predictions(validation, :);
        %classifications.confidences(validation, :, :) = tmp.confidences(validation, :, :);
        classifications.correctness(validation, :) = tmp.correctness(validation, :);
        classifications.correct_total = classifications.correct_total + tmp.correctness(validation, :);
    end
    
    % Summarise across validations
    classifications.accuracy_per_class = classifications.correct_total / size(classifications.leave_outs, 1);
    classifications.accuracy = sum(classifications.correct_total) / numel(classifications.leave_outs);
    
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
            tmp = svm_lol_libsvm(values, validation);
            
            % Store values used and results
            classifications{fly}.class_data(validation, :) = tmp.class_data;
            classifications{fly}.svms{validation} = tmp.svms{validation};
            classifications{fly}.predictions(validation, :) = tmp.predictions(validation, :);
            %classifications{fly}.confidences(validation, :, :) = tmp.confidences(validation, :, :);
            classifications{fly}.correctness(validation, :) = tmp.correctness(validation, :);
            classifications{fly}.correct_total = classifications{fly}.correct_total + tmp.correctness(validation, :);
        end
        
        % Summarise across validations
        classifications{fly}.accuracy_per_class = classifications{fly}.correct_total / size(classifications{fly}.leave_outs, 1);
        classifications{fly}.accuracy = sum(classifications{fly}.correct_total) / numel(classifications{fly}.leave_outs);
    end
    
end

%% Save

%save([results_location results_file], 'classifications');

%%
% 
% figure;
% 
% foo = cat(1, values(:, :, 1), values(:, :, 2));
% [coeff, score, latent, ~, explained] = pca(foo);
% 
% scatter(score((1:13), 1), score((1:13), 2), 'r'); hold on;
% scatter(score((1:13)+13, 1), score((1:13)+13, 2), 'b');
% 
% %%
% 
% figure;
% 
% for fly = 1 : size(values_all, 3)
%     subplot(size(values_all, 3), 1, fly);
%     
%     bar(squeeze(values_all(:, :, fly, 1) - values_all(:, :, fly, 2)));
% end


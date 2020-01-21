function [classification] = svm_loo_liblinear_manual(values, cost_range)
% Classifies using SVM (liblinear)
%
% Conducts nested cross-validation with leave-one-out at both stages
%   Outer stage - repeated for each left-out observation
%   Inner stage - repeated for each left-out observation, out of the
%       remaining observations
%   Inner stage - repeated at each cost-value, validated on left-out
%       observation
%   Outer stage - tests left-out observation on best cost-value model from
%       inner stage
%
% Normalises data 
%
% Assumes equal number of observations in each class
%
% Inputs:
%   values: 3D matrix (observations x features x classes)
%   cost_range: vector; list of costs to try; actual costs, e.g. 2^-10 (not a power to apply, in this case -10)
% Outputs:
%   classification: struct holding classification results

nObservations = size(values, 1);
nFeatures = size(values, 2);
nClasses = size(values, 3);

% % Normalise all data
% values_normed = zeros(nObservations * nClasses, nFeatures);
% observation_counter = (1:nObservations);
% for class = 1 : nClasses
%     values_normed(observation_counter, :) = values(:, :, class);
%     observation_counter = observation_counter + nObservations;
% end
% values_normed = zscore(values_normed, [], 1);
% observation_counter = (1:nObservations);
% for class = 1 : nClasses
%     values(:, :, class) = values_normed(observation_counter, :);
%     observation_counter = observation_counter + nObservations;
% end

% Build feature-observation matrix for each class, and labels
class_data = cell(1, nClasses);
class_labels = cell(1, nClasses);
for class = 1 : nClasses
    class_data{class} = values(:, :, class);
    class_labels{class} = zeros(nObservations, 1) + class;
end

% Leave-one-out combinations
% Leave out first of each class, then second, then third, etc
% Leave-outs x classes
% Each row gives combo of observations (which observations from which class to leave out
leave_outs = repmat((1:nObservations)', [1 nClasses]);

% Training and test labels will always be the same
data_train_labels = zeros((nObservations-1) * nClasses, 1);
data_test_labels = zeros(nClasses, 1);
observation_counter = (1:nObservations-1);
for class = 1 : nClasses
    data_train_labels(observation_counter) = class;
    data_test_labels(class) = class;
    observation_counter = observation_counter + (nObservations-1);
end

% Training and test labels will always be the same (validation)
data_vtrain_labels = zeros((nObservations-2) * nClasses, 1);
data_vtest_labels = zeros(nClasses, 1);
observation_counter = (1:nObservations-2);
for class = 1 : nClasses
    data_vtrain_labels(observation_counter) = class;
    data_vtest_labels(class) = class;
    observation_counter = observation_counter + (nObservations-2);
end

% LIBLINEAR requires double for labels
data_train_labels = double(data_train_labels);
data_test_labels = double(data_test_labels);
data_vtrain_labels = double(data_vtrain_labels);
data_vtest_labels = double(data_vtest_labels);

% Counter of classification correctness
correct_total = zeros(1, nClasses);

% Storage of SVMs and classifications
svms = cell(size(leave_outs, 1), 1);
predictions = zeros(size(leave_outs));
validation_costs = zeros(size(leave_outs, 1), 1);
validation_cost_ids = zeros(size(leave_outs, 1), 1);
correctness = zeros(size(leave_outs));

% Default cost = 2^0
if nargin == 1
    cost_range = 1;
end

% Loop through each leave_out combination
for test_point = 1 : size(values, 1)
    
    values_validate = values;
    values_validate(test_point, :, :) = [];
    
    cost_vAccuracies = zeros(length(cost_range), 1);
    
    % Find best hyperparameter value (cost) using cross-validation
    for cost_counter = 1 : length(cost_range)
        cost = cost_range(cost_counter);
        
        vAccuracies = zeros(size(values_validate, 1), 1);
        
        for validate_point = 1 : size(values_validate, 1)
            
            % Get training and validation data
            data_vtrain = values_validate;
            data_vtrain(validate_point, :, :) = [];
            data_vtest = values_validate(validate_point, :, :);
            
            % Format into (observations*classes x features)
            data_vtrain = reshape(permute(data_vtrain, [1 3 2]), [size(data_vtrain, 1)*size(data_vtrain, 3) size(data_vtrain, 2)]);
            data_vtest = reshape(permute(data_vtest, [1 3 2]), [size(data_vtest, 1)*size(data_vtest, 3) size(data_vtest, 2)]);
            
            % Normalise training data
            means = mean(data_vtrain, 1);
            stds = std(data_vtrain, [], 1);
            means_mat = repmat(means, [size(data_vtrain, 1), 1]);
            stds_mat = repmat(stds, [size(data_vtrain, 1), 1]);
            data_vtrain = (data_vtrain - means_mat) ./ stds_mat;
            
            % Normalise validation data (using same parameters as for
            % training data)
            means_mat = repmat(means, [size(data_vtest, 1), 1]);
            stds_mat = repmat(stds, [size(data_vtest, 1), 1]);
            data_vtest = (data_vtest - means_mat) ./ stds_mat;
            
            %data_vtrain(isnan(data_vtrain)) = 0;
            %data_vtest(isnan(data_vtest)) = 0;
            
            % Train SVM
            trained = train(data_vtrain_labels, sparse(data_vtrain), ['-c ' num2str(cost) ' -q']);
            
            % Classify validation data using SVM
            [prediction, accuracy, confidence] = predict(data_vtest_labels, sparse(data_vtest), trained, '-q');
            
            vAccuracies(validate_point) = accuracy(1);
            
        end
        
        cost_vAccuracies(cost_counter) = mean(vAccuracies);
        
    end
    
    % Identify best cost from validation set
    [validation_accuracy, cost_counter] = max(cost_vAccuracies);
    validation_cost = cost_range(cost_counter);
    
    % Training data is the validation data
    data_train = values_validate;
    
    % Testing data
    data_test = values(test_point, :, :);
    
    % Format into (observations*classes x features)
    data_train = reshape(permute(data_train, [1 3 2]), [size(data_train, 1)*size(data_train, 3) size(data_train, 2)]);
    data_test = reshape(permute(data_test, [1 3 2]), [size(data_test, 1)*size(data_test, 3) size(data_test, 2)]);
    
    % Normalise all validation data
    means = mean(data_train, 1);
    stds = std(data_train, [], 1);
    means_mat = repmat(means, [size(data_train, 1), 1]);
    stds_mat = repmat(stds, [size(data_train, 1), 1]);
    data_train = (data_train - means_mat) ./ stds_mat;
    
    % Normalise test data (using same parameters as for validation data)
    means_mat = repmat(means, [size(data_test, 1), 1]);
    stds_mat = repmat(stds, [size(data_test, 1), 1]);
    data_test = (data_test - means_mat) ./ stds_mat;
    
%     data_train(isnan(data_train)) = 0;
%     data_test(isnan(data_test)) = 0;
    
    % Train SVM with cost identified from validation data
    trained = train(data_train_labels, sparse(data_train), ['-c ' num2str(validation_cost) ' -q']);
    
    % Classify training data using SVM
    [prediction, accuracy, confidence] = predict(data_test_labels, sparse(data_test), trained, '-q');
    
    % Store
    svms{test_point} = trained;
    predictions(test_point, :) = prediction;
    validation_accuracies(test_point) = validation_accuracy;
    validation_costs(test_point) = validation_cost;
    validation_cost_ids(test_point) = cost_counter;
    
    % Check accuracy of classifications
    correct = prediction == data_test_labels;
    correctness(test_point, :) = correct;
    correct_total = correct_total + correct';
    
end

accuracy_per_class = correct_total / size(leave_outs, 1);
accuracy = sum(correct_total) / numel(leave_outs);

classification.class_data = class_data;
classification.leave_outs = leave_outs;
classification.svms = svms;
classification.costs = cost_range;
classification.validation_accuracies = validation_accuracies;
classification.validation_costs = validation_costs;
classification.validation_cost_ids = validation_cost_ids;
classification.predictions = predictions;
%classification.confidences = confidences;
classification.correctness = correctness;
classification.correct_total = correct_total;
classification.accuracy_per_class = accuracy_per_class;
classification.accuracy = accuracy;

end


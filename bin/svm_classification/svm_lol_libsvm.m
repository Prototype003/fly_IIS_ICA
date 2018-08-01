function [classification] = svm_lol_libsvm(values, leave_out_counter)
% Classifies using SVM, validating using leave-one-out
% Assumes equal number of observations in each class
%
% Inputs:
%   values: 3D matrix (observations x features x classes)
%   leave_out_counter: optional int; counter for which data-combo to leave out
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

% Generate comprehensive leave-one-out combinations
observation_labels = repmat((1:nObservations), [1 nClasses]);
leave_outs = unique(nchoosek(observation_labels, nClasses), 'rows'); % Each row gives a unique combo of observations to leave out

% Training and test labels will always be the same
data_train_labels = zeros((nObservations-1) * nClasses, 1);
data_test_labels = zeros(nClasses, 1);
observation_counter = (1:nObservations-1);
for class = 1 : nClasses
    data_train_labels(observation_counter) = class;
    data_test_labels(class) = class;
    observation_counter = observation_counter + (nObservations-1);
end

% Counter of classification correctness
correct_total = zeros(1, nClasses);

% Storage of SVMs and classifications
svms = cell(size(leave_outs, 1), 1);
predictions = zeros(size(leave_outs));
confidences = zeros([size(leave_outs) nClasses]);
correctness = zeros(size(leave_outs));

if nargin > 1
    % Perform the specified leave out training/testing
    
    leave_out = leave_outs(leave_out_counter, :);
    
    % Get data, with leave_out observations removed
    data_train = zeros((nObservations-1)*nClasses, nFeatures);
    data_test = zeros(length(leave_out), nFeatures);
    observation_counter = (1:nObservations-1);
    for class = 1 : nClasses
        data_train(observation_counter, :) = class_data{class}(1:end ~= leave_out(class), :); % training data
        data_test(class, :) = class_data{class}(leave_out(class), :); % test data
        observation_counter = observation_counter + (nObservations-1);
    end
    
%     % Normalise training data ([0 1])
%     % https://stackoverflow.com/questions/17335738/results-on-libsvm-favors-only-one-class-from-two-classes
%     % https://stats.stackexchange.com/questions/77350/perform-feature-normalization-before-or-within-model-validation
%     % https://stats.stackexchange.com/questions/70801/how-to-normalize-data-to-0-1-range
%     max_values = max(data_train, [], 1); % Max across rows (per column)
%     min_values = min(data_train, [], 1); % Min across rows (per column)
%     max_mat = repmat(max_values, [size(data_train, 1) 1]);
%     min_mat = repmat(min_values, [size(data_train, 1) 1]);
%     data_train = (data_train - min_mat) ./ (max_mat - min_mat);
%     
% %     % Normalise testing data (using same parameters as for training data)
% %     max_mat = repmat(max_values, [size(data_test, 1) 1]);
% %     min_mat = repmat(min_values, [size(data_test, 1) 1]);
% %     data_test = (data_test - min_mat) ./ (max_mat - min_mat);
%     
%     % Normalise testing data independently
%     max_values = max(data_test, [], 1); % Max across rows (per column)
%     min_values = min(data_test, [], 1); % Min across rows (per column)
%     max_mat = repmat(max_values, [size(data_test, 1) 1]);
%     min_mat = repmat(min_values, [size(data_test, 1) 1]);
%     data_test = (data_test - min_mat) ./ (max_mat - min_mat);
    
    % Normalise training data ([-1 1])
    means = mean(data_train, 1);
    stds = std(data_train, [], 1);
    means_mat = repmat(means, [size(data_train, 1), 1]);
    stds_mat = repmat(stds, [size(data_train, 1), 1]);
    data_train = (data_train - means_mat) ./ stds_mat;
    
    % Normalise testing data (using same parameters as for training data)
    means_mat = repmat(means, [size(data_test, 1), 1]);
    stds_mat = repmat(stds, [size(data_test, 1), 1]);
    data_test = (data_test - means_mat) ./ stds_mat;
    
    % Normalise data (training and testing data independently)
    %data_train = zscore(data_train, [], 1);
    %data_test = zscore(data_test, [], 1);
    
    % Train SVM
    trained = svmtrain(data_train_labels, data_train, '-t 0'); % '-t 0' = linear; '-t 2' = rbf ('-t 2' is default)
    
    % Classify test data using SVM
    [prediction, accuracy, confidence] = svmpredict(data_test_labels, data_test, trained);
    
    % Store
    svms{leave_out_counter} = trained;
    predictions(leave_out_counter, :) = prediction;
    %confidences(leave_out_counter, :, :) = confidence;
    
    % Check accuracy of classifications
    correct = prediction == data_test_labels;
    correctness(leave_out_counter, :) = correct;
    correct_total = correct_total + correct';
    
else
    % Loop through each leave_out combination
    for leave_out_counter = 1 : size(leave_outs, 1)
        
        leave_out = leave_outs(leave_out_counter, :);
        
        % Get data, with leave_out observations removed
        data_train = zeros((nObservations-1)*nClasses, nFeatures);
        data_test = zeros(length(leave_out), nFeatures);
        observation_counter = (1:nObservations-1);
        for class = 1 : nClasses
            data_train(observation_counter, :) = class_data{class}(1:end ~= leave_out(class), :); % training data
            data_test(class, :) = class_data{class}(leave_out(class), :); % test data
            observation_counter = observation_counter + (nObservations-1);
        end
        
        % Normalise training data ([-1 1])
        means = mean(data_train, 1);
        stds = std(data_train, [], 1);
        means_mat = repmat(means, [size(data_train, 1), 1]);
        stds_mat = repmat(stds, [size(data_train, 1), 1]);
        data_train = (data_train - means_mat) ./ stds_mat;
        
        % Normalise testing data (using same parameters as for training data)
        means_mat = repmat(means, [size(data_test, 1), 1]);
        stds_mat = repmat(stds, [size(data_test, 1), 1]);
        data_test = (data_test - means_mat) ./ stds_mat;
        
        % Train SVM
        trained = svmtrain(data_train_labels, data_train, '-t 0'); % '-t 0' = linear; '-t 2' = rbf ('-t 2' is default)
        
        % Classify test data using SVM
        [prediction, accuracy, confidence] = svmpredict(data_test_labels, data_test, trained);
        
        % Store
        svms{leave_out_counter} = trained;
        predictions(leave_out_counter, :) = prediction;
        %confidences(leave_out_counter, :, :) = confidence;
        
        % Check accuracy of classifications
        correct = prediction == data_test_labels;
        correctness(leave_out_counter, :) = correct;
        correct_total = correct_total + correct';
        
    end
end

accuracy_per_class = correct_total / size(leave_outs, 1);
accuracy = sum(correct_total) / numel(leave_outs);

classification.class_data = class_data;
classification.leave_outs = leave_outs;
classification.svms = svms;
classification.predictions = predictions;
%classification.confidences = confidences;
classification.correctness = correctness;
classification.correct_total = correct_total;
classification.accuracy_per_class = accuracy_per_class;
classification.accuracy = accuracy;

end


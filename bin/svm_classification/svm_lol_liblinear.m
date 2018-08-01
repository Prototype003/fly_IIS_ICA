function [accuracy] = svm_lol_liblinear(values)
% Classifies using SVM, validating using leave-one-out
% Assumes equal number of observations in each class
%
% Inputs:
%   values: 3D matrix (observations x features x classes)
%   leave_out_counter: optional int; counter for which data-combo to leave out
% Outputs:
%   accuracy: accuracy after searching for cost parameter C and k-fold validation

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
data_train_labels = zeros((nObservations) * nClasses, 1);
observation_counter = (1:nObservations);
for class = 1 : nClasses
    data_train_labels(observation_counter) = class;
    observation_counter = observation_counter + (nObservations);
end

% Get data, with leave_out observations removed
data_train = zeros((nObservations)*nClasses, nFeatures);
observation_counter = (1:nObservations);
for class = 1 : nClasses
    data_train(observation_counter, :) = class_data{class}(1:end, :); % training data
    observation_counter = observation_counter + (nObservations);
end

c = train(data_train_labels, sparse(data_train), ['-C']);
c_param = c(1);

% Train SVM
result = train(data_train_labels, sparse(data_train), ['-c ' num2str(c_param) ' -v 16']); % Switching from -s 2 to -s 1 is fine

accuracy = result(1);

% % Store
% svms{leave_out_counter} = trained;
% predictions(leave_out_counter, :) = prediction;
% 
% 
% accuracy_per_class = correct_total / size(leave_outs, 1);
% accuracy = sum(correct_total) / numel(leave_outs);
% 
% classification.class_data = class_data;
% classification.leave_outs = leave_outs;
% classification.svms = svms;
% classification.predictions = predictions;
% %classification.confidences = confidences;
% classification.correctness = correctness;
% classification.correct_total = correct_total;
% classification.accuracy_per_class = accuracy_per_class;
% classification.accuracy = accuracy;

end


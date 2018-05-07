function [accuracy, leave_out, classifications, correctness] = classify_lol(class_data)
% Trains a 1D classifier using the leave-one-out method
%
% Inputs:
%   class_data: matrix (data points x classes)
%
% Outputs:
%   accuracy: number; average accuracy across all leave-out combinations
%   leave_out: matrix (datapoint combinations between classes x classes); each row corresponds to a datapoint pair (one point from each class) which
%       is left out from training for classification (number of classifications conducted = number of rows in this matrix)
%   classifications: matrix (same size as leave_out); holds the class assigned to the left out datapoints
%   correctness: matrix (same size as leave_out); holds the correctness of classification to the left out datapoints; 1 = correct, 0 = incorrect
%

[nPoints, nClasses] = size(class_data);

% Refer to these when getting data vectors with a point left out
class_data_vectors = cell(1, nClasses);
for class = 1 : nClasses
    class_data_vectors{class} = class_data(:, class);
end

% Combinations of class datapoints to leave out (e.g. [x1,y1], [x1,y2], ... [xn,yn])
leave_out = combvec((1:nPoints), (1:nPoints))'; % TODO: update so it is dynamic for nClasses

% Classification results storage
classifications = zeros(size(leave_out, 1), nClasses);
correctness = zeros(size(leave_out, 1), nClasses);

% Iterate through combinations of datapoints from each class (e.g. [x1,y1], [x1,y2], ... [xn,yn]) and train, classify
class_data_lol = zeros(nPoints-1, nClasses);
for exclusion = 1 : size(leave_out, 1)
    
    % Exclude datapoints
    for class = 1 : nClasses
        include_points = true(nPoints, 1);
        include_points(leave_out(exclusion, class)) = 0;
        class_data_lol(:, class) = class_data_vectors{class}(include_points);
    end
    
    % Find class centers after exclusion
    class_centers = mean(class_data_lol, 1);
    
    % Determine which class each excluded datapoint is closest to
    for class = 1 : nClasses % class is the correct class
        class_distances = abs(class_centers - class_data(leave_out(exclusion, class), class));
        [~, classifications(exclusion, class)] = min(class_distances);
        correctness(exclusion, class) = classifications(exclusion, class) == class;
    end
end

accuracy = 100 * sum(correctness(:)) / numel(correctness);

end


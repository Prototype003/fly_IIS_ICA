function [deviation_total] = error_distance(data, type)
% Calculates mean of each class, and deviation of each 'error' datapoint from their correct mean
% Only coded for 2 classes
%
% Inputs:
%   data: n x classes matrix; each column holds the datapoints for a class;
%       classes should be ordered so the first column has in principle the highest values, second column the second highest, etc
%   type: string; 'para' for parametric deviation or 'nonpara' for non-parametric
%
% Outputs:
%   deviation: sum of error deviations

% Label datapoints with their class
data_classes = ones(size(data));
for class = 2 : size(data, 2)
    data_classes(:, class) = class;
end

if strcmp(type, 'para')
    deviation_total = para(data, data_classes);
elseif strcmp(type, 'nonpara')
    deviation_total = nonpara(data);
else
    % Do nothing x.x
end

    function [deviation_total] = para(data, data_classes)
        % Parametric deviation
        
        % Order datapoints across classes
        [data_sorted, sort_order] = sort(data(:), 'desc');
        data_classes = data_classes(sort_order);
        
        % Separate into 'new' posterior groups; each column is a class
        data_sorted = reshape(data_sorted, size(data));
        data_classes = reshape(data_classes, size(data));
        
        % Get distributions based on ordered values
        data_sorted_means = mean(data_sorted, 1);
        data_sorted_stds = std(data_sorted, [], 1);
        pooled_std = sqrt(sum(data_sorted_stds.^2));
        
        % Find errors in the new groups and sum their deviations to their correct group
        deviation_total = 0;
        for class = 1 : size(data, 2)
            for point = 1 : size(data, 1)
                if class ~= data_classes(point, class)
                    % the datapoint is in the wrong class
                    deviation_total = deviation_total + (abs(data_sorted(point, class) - data_sorted_means(class)) / pooled_std);
                end
            end
        end
        
    end

    function [deviation_total] = nonpara(data)
        % Non-parametric deviation (rank correlation with ideal ordering)
        
        % Order within class
        class_sorted = zeros(size(data));
        for class = 1 : size(data, 2)
            class_sorted(:, class) = sort(data(:, class), 'desc');
        end
        ideal_ranking = reshape((1:numel(data)), size(data));
        
        % Order across classes
        [data_sorted, actual_ranking] = sort(class_sorted(:), 'desc');
        
        deviation_total = corr(ideal_ranking(:), actual_ranking, 'Type', 'Kendall'); % 'Kendall' or 'Spearman' for rank order correlation
        deviation_total = corr(class_sorted(:), data_sorted, 'Type', 'Kendall'); % Should be the same as when using rankings (because the correlation turns values into rankings anyways)
        
    end

end


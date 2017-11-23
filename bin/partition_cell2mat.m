function [partition_matrix] = partition_cell2mat(partition_cell)
% Converts a cell partition (partition{1}--/-->partition{2} to a matrix
% (partition(1, :)--/-->partition(2, :)). Uneven partitions are made even
% with -1 (e.g. (1,)--/-->(2,3) = (1,-1)--/-->(2,3))
%
% Inputs:
%   partition_cell - 
%
% Outputs:
%   partition_matrix - 

from = partition_cell{1};
to = partition_cell{2};

partition_size = max([length(from), length(to)]);
partition_matrix = zeros(2, partition_size);
partition_matrix = partition_matrix - 1; % Fill uneven groups with placeholder -1

partition_matrix(1, 1:length(from)) = from;
partition_matrix(2, 1:length(to)) = to;

end


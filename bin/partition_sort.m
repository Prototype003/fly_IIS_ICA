function [sort_indices] = partition_sort(partitions)
% Hardcoded sort

% Desired partition order
partition_order = {...
    {[0], [1 2]},...
    {[1 2], [0]},...
    {[0 1], [2]},...
    {[2], [0 1]},...
    {[0 2], [1]},...
    {[1], [0 2]},...
    };

sort_indices = zeros(1, length(partition_order));

for partition_counter = 1 : length(partition_order)
    partition = partition_order{partition_counter};
    for partitions_counter = 1 : length(partitions)
        if isequal(partitions{partitions_counter}, partition)
            sort_indices(partition_counter) = partitions_counter;
            break
        end
    end
end

end
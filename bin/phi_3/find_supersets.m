%% Function: find channel sets which contain some subset
function [superset_ids] = find_supersets(subset, supersets)
% Returns indexes of sets which contain some given subset
%
% Inputs:
%   subset: 1D vector;
%   supersets: 2D matrix; each row gives a superset; size of second dimension
%       should be larger than the size of 'subset'
%
% Outputs:
%   superset_ids: row indices of supersets which contain the subset

superset_status = zeros(size(supersets, 1), 1);
for superset = 1 : size(supersets, 1)
    superset_status(superset) = all(ismember(subset, supersets(superset, :)));
end

superset_ids = find(superset_status);

end
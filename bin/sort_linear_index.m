function [sortedA,sortIndex] = sort_linear_index(A,sortDim,sortOrder)
%#SORT_LINEAR_INDEX   Just like SORT, but returns linear indices

  sizeA = size(A);  %# Get the matrix size
  if nargin < 2
    sortDim = find(sizeA > 1,1);  %# Define sortDim, if necessary
  end
  if nargin < 3
    sortOrder = 'ascend';  %# Define sortOrder, if necessary
  end
  [sortedA,sortIndex] = sort(A,sortDim,sortOrder);  %# Sort the matrix
  [subIndex{1:numel(sizeA)}] = ...  %# Create a set of matrix subscripts
     ind2sub(sizeA,reshape(1:prod(sizeA),sizeA));
  subIndex{sortDim} = sortIndex;  %# Overwrite part of the subscripts with
                                  %#   the sort indices
  sortIndex = sub2ind(sizeA,subIndex{:});  %# Find the linear indices

end
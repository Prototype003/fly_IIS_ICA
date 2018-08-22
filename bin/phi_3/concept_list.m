function [concept_list_full] = concept_list(nChannels)
% Returns full list of possible concepts (powerset of 'nChannels' elements)
%
% Inputs:
%   nChannels: int, how many channels to consider (elements in the system)
% Outputs:
%   concept_list_full: cell array; each cell holds a vector representing each
%       possible concept; 1st order concepts come first, sorted by element,
%       then 2nd order concepts, 3rd order concepts, etc.
%       NOTE: elements are 1-indexed (remember when joining results from
%       python that elements are 0-indexed)

nConcepts = 0; % power-set of nChannels

concept_list_full = cell(0); concept_list_counter = 1;
for concept_order = 1 : nChannels
    subsets = nchoosek((1:nChannels), concept_order);
    nConcepts = nConcepts + size(subsets, 1);
    for concept = 1 : size(subsets, 1)
        concept_list_full{concept_list_counter} = subsets(concept, :); % - 1; % '-1' if 0-indexed
        concept_list_counter = concept_list_counter + 1;
    end
end

end
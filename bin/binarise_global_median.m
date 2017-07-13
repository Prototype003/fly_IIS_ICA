%% Function: global median split

function [binarised] = binarise_global_median(fly_data)
% MATLAB conversion of binarise_global_median() in fly_phi.py
% Finds the median for each channel, across all epoch-trials
% Binarises samples based on the median: 1 if greater than median, 0 otherwise

% Reorder dimensions to prepare for reshape
reordered = permute(fly_data, [1 3 2 4 5]);
% Reformat - collapse on trials so that the first dimension is samples x trials
sizes = size(reordered);
globalised = reshape(reordered, [sizes(1)*sizes(2) sizes(3) sizes(4) sizes(5)]);

% Get median per channel, fly, condition
channel_medians = median(globalised, 1);

% Binarise based on median:
% Greater than median = 1
% Less than or equal to median = 0
channel_medians = repmat(channel_medians, [sizes(1)*sizes(2) 1 1 1]); % Repeat median for each sample, for comparison
binarised = globalised > channel_medians;

% Return to original shape
binarised = reshape(binarised, [sizes(1) sizes(2) sizes(3) sizes(4) sizes(5)]);
binarised = permute(binarised, [1 3 2 4 5]);

end
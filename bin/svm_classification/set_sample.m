function [features] = set_sample(channel_sets, nFeatures)
% Randomly selects n sets, such that across n sets all channels are used at least once
%
% Inputs:
%   channel_sets: matrix of potential sets (sets x channels)
%   nFeatures: integer, number of sets to pick out
%
% Outputs:
%   features: vector of indices to the selected sets (length = nFeatures)

nChannels = max(channel_sets);

if nFeatures * size(channel_sets, 2) < nChannels
    disp('Not enough features to cover all channels');
    return
end

% Keep original copy so we can reset if required conditions aren't met
sample_space = channel_sets;

% Get random sample
features = randperm(size(channel_sets, 1), nFeatures);

% Check if selected sets cover all channels
covered = 0;
channels = channel_sets(features, :);
channels = unique(channels(:));
if numel(channels) == max(nChannels)
    covered = 1;
end

% Repeat while not all channels are covered
while covered == 0
    
    % Keep original copy so we can reset if required conditions aren't met
    sample_space = channel_sets;
    
    % Get random sample
    features = randperm(size(channel_sets, 1), nFeatures);
    
    % Check if selected sets cover all channels
    covered = 0;
    channels = channel_sets(features, :);
    channels = unique(channels(:));
    if numel(channels) == max(nChannels)
        covered = 1;
    end
    
end

end


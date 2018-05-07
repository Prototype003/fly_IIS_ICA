function [phi_star, H, H_cond, MI, MI_star, MIP, partitions, partitions_phis, partitions_H, partitions_H_cond, partitions_MI, partitions_MI_star] = phistar_mip(cov_past_past, cov_past_present, cov_present_present, channels)
% Finds phi-star and associated metrics corresponding to the MIP
% Uses Hauns' Feb2014 toolbox
%
% To do this, calculates phi-star for every possible partition
%
% Inputs:
%   Channels = vector of channels (node labels)
%
% Outputs:

nChannels = length(channels);
channel_map = (1:nChannels);

% Get all possible partitions
partitions = SetPartition(channel_map);
partitions = partitions(2:end); % The first partitioning is the unpartitioned system, so we ignore that

% Setup storage of values for each partition
partitions_phis = zeros(size(partitions, 1), 2);
partitions_H = zeros(size(partitions));
partitions_H_cond = zeros(size(partitions));
partitions_MI = zeros(size(partitions));
partitions_MI_star = zeros(size(partitions));

% For each partition, find phi-star
% MIP is where phi-star is the smallest

% Get phi for the first partition as baseline
partition = partitions{1};
partition_formatted = zeros(nChannels, 1);
for group = 1 : length(partition)
    partition_formatted(partition{group}) = group;
end
[phi_star, ~, ~, ~, MI, ~, ~, H_cond, H, MI_star] = phi_comp(cov_present_present, cov_past_present', cov_past_past, [], partition_formatted);
MIP = partition;
% Save values for partition
partitions_phis(1, :) = phi_star;
partitions_H(1) = H;
partitions_H_cond(1) = H_cond;
partitions_MI(1) = MI;
partitions_MI_star(1) = MI_star;

% Repeat for each partition and update if phi is smaller
for partition_counter = 2 : length(partitions) % Already checked the first partition, so start from the second
    partition = partitions{partition_counter};
    
    % Reformat partition for input to phi_star_Gauss()
    partition_formatted = zeros(nChannels, 1);
    for group = 1 : length(partition)
        partition_formatted(partition{group}) = group;
    end
    
    [phi_star_new, ~, ~, ~, MI_new, ~, ~, H_cond_new, H_new, MI_star_new] = phi_comp(cov_present_present, cov_past_present', cov_past_past, [], partition_formatted);
    
    if phi_star_new(2) < phi_star(2)
        MIP = partition;
        phi_star = phi_star_new;
        H = H_new;
        H_cond = H_cond_new;
        MI = MI_new;
        MI_star = MI_star_new;
    end
    
    % Save values for partition
    partitions_phis(partition_counter, :) = phi_star_new;
    partitions_H(partition_counter) = H_new;
    partitions_H_cond(partition_counter) = H_cond_new;
    partitions_MI(partition_counter) = MI_new;
    partitions_MI_star(partition_counter) = MI_star_new;
end

% MIP is mapped (to consecutive channels), so unmap (to original labels)
for group = 1 : length(MIP)
    for channel_counter = 1 : length(MIP{group})
        channel = channels(MIP{group}(channel_counter));
        MIP{group}(channel_counter) = channel;
    end
end


end


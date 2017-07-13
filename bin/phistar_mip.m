function [phi_star, H, H_cond, MI, MI_star, MIP] = phistar_mip(cov_past_past, cov_past_present, cov_present_present, channels)
% Finds phi-star and associated metrics corresponding to the MIP
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

% For each partition, find phi-star
% MIP is where phi-star is the smallest

% Get phi for the first partition as baseline
partition = partitions{2}; % The first partitioning is the unpartitioned system, so we ignore that
partition_formatted = zeros(nChannels, 1);
for group = 1 : length(partition)
    partition_formatted(partition{group}) = group;
end
[phi_star, ~, ~, ~, MI, ~, ~, H_cond, H, MI_star] = phi_comp(cov_present_present, cov_past_present, cov_past_past, [], partition_formatted);
MIP = partition;

% Repeat for each partition and update if phi is smaller
for partition_counter = 3 : length(partitions) % We're ignoring the 1st partition (unpartitioned), and we've already done the 2nd partitioning, so start from the 3rd
    partition = partitions{partition_counter};
    
    % Reformat partition for input to phi_star_Gauss()
    partition_formatted = zeros(nChannels, 1);
    for group = 1 : length(partition)
        partition_formatted(partition{group}) = group;
    end
    
    [phi_star_new, ~, ~, ~, MI_new, ~, ~, H_cond_new, H_new, MI_star_new] = phi_comp(cov_present_present, cov_past_present, cov_past_past, [], partition_formatted);
    
    if phi_star_new(2) < phi_star(2)
        MIP = partition;
        phi_star = phi_star_new;
        H = H_new;
        H_cond = H_cond_new;
        MI = MI_new;
        MI_star = MI_star_new;
    end
end

% MIP is mapped (to consecutive channels), so unmap (to original labels)
for group = 1 : length(MIP)
    for channel_counter = 1 : length(MIP{group})
        channel = channels(MIP{group}(channel_counter));
        MIP{group}(channel_counter) = channel;
    end
end


end


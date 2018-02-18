function [phi_g, MIP, partitions, partitions_phis] = phig_mip(cov_past_past, cov_past_present, cov_present_present, channels)
% Finds phi-g and associated metrics corresponding to the MIP
%
% To do this, calculates phi-g for every possible partition
%
% Based off phistar_mip.m
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

% For each partition, find phi-g
% MIP is where phi-g is the smallest

% Get phi for the first partition as baseline
partition = partitions{1};
partition_formatted = zeros(nChannels, 1);
for group = 1 : length(partition)
    partition_formatted(partition{group}) = group;
end
%[phi_star, ~, ~, ~, MI, ~, ~, H_cond, H, MI_star] = phi_comp(cov_present_present, cov_past_present, cov_past_past, [], partition_formatted);
[phi_g] = phi_G_Gauss(cov_past_past, cov_past_present, cov_present_present, partition_formatted); phi_g = [phi_g phi_g/phi2_normFactor(partition_formatted, cov_past_past)];
MIP = partition;
% Save values for partition
partitions_phis(1, :) = phi_g;

% Repeat for each partition and update if phi is smaller
for partition_counter = 2 : length(partitions) % Already checked the first partition, so start from the second
    partition = partitions{partition_counter};
    
    % Reformat partition for input to phi_star_Gauss()
    partition_formatted = zeros(nChannels, 1);
    for group = 1 : length(partition)
        partition_formatted(partition{group}) = group;
    end
    
    %[phi_star_new, ~, ~, ~, MI_new, ~, ~, H_cond_new, H_new, MI_star_new] = phi_comp(cov_present_present, cov_past_present, cov_past_past, [], partition_formatted);
    [phi_g_new] = phi_G_Gauss(cov_past_past, cov_past_present, cov_present_present, partition_formatted); phi_g_new = [phi_g_new phi_g_new/phi2_normFactor(partition_formatted, cov_past_past)];
    
    if phi_g_new(2) < phi_g(2)
        MIP = partition;
        phi_g = phi_g_new;
    end
    
    % Save values for partition
    partitions_phis(partition_counter, :) = phi_g_new;
end

% MIP is mapped (to consecutive channels), so unmap (to original labels)
for group = 1 : length(MIP)
    for channel_counter = 1 : length(MIP{group})
        channel = channels(MIP{group}(channel_counter));
        MIP{group}(channel_counter) = channel;
    end
end


end


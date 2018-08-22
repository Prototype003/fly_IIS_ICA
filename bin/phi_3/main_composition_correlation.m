%% Description

%{

For seeing if there's a relationship between two-channel phi and their
"corresponding" mechanism within 4-channel phi

%}

%% Setup

nChannels = [2 4];

%% Load relevant files

% global_tpm = 0;
% 
% source_dir = 'results/';
% 
% phis = cell(4, 1);
% 
% % 2 channels
% source_2ch = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels2_globalTPM' num2str(global_tpm) '.mat'];
% tmp = load([source_dir source_2ch]);
% phis{2} = tmp.phis{1};
% phis{2}.concept_list_full = concept_list(2); % Relevant 2nd order concept is the 3rd one
% 
% % 4 channels
% tic
% disp('loading 4ch');
% source_4ch = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM' num2str(global_tpm) '.mat'];
% tmp = load([source_dir source_4ch]);
% phis{4} = tmp.phis{1};
% disp('loaded 4ch');
% phis{4}.concept_list_full = concept_list(4);
% toc

%% Compare second-order concepts

% For each 2-channel set, find all 4-channel sets which contain the two channels
% Then get the 2nd order small-phi corresponding to the 2 channels

% Find how many data entries we need
matches = 0;
for channel_pair_counter = 1 : size(phis{2}.channel_sets, 1)
    channel_pair = phis{2}.channel_sets(channel_pair_counter, :);
    
    supersets = find_supersets(channel_pair, phis{4}.channel_sets);
    matches = matches + length(supersets);
end

% Storage structure
tmp = size(phis{4}.big_mips);
storage_size = [matches tmp([5 6 7 8 2])];
concept_phi = struct();
concept_phi.twoCh = zeros(storage_size);
concept_phi.fourCh = zeros(storage_size);

match_counter = 1; 
for channel_pair_counter = 1 : size(phis{2}.channel_sets, 1)
    channel_pair = phis{2}.channel_sets(channel_pair_counter, :);
    
    supersets = find_supersets(channel_pair, phis{4}.channel_sets);
    
    for superset_counter = 1 : length(supersets)
        superset = supersets(superset_counter);
        
        % Determine which concept should be compared
        %   Keep in mind that which 2nd order concept to use will depend on
        %   the position of the channel-pair among the 4 channels
        
        % Convert channel IDs to IDs within the 4-channel set
        channel_pair_location = zeros(size(channel_pair));
        for channel_counter = 1 : length(channel_pair)
            channel = channel_pair(channel_counter);
            channel_pair_location(channel_counter) = find(phis{4}.channel_sets(superset, :) == channel);
        end
        
        % Get small-phi corresponding to the relevant concept
        for concept = 1 : length(phis{4}.concept_list_full)
            if isequal(phis{4}.concept_list_full{concept}, channel_pair_location)
                % 2ch - get weighted average of small phi across states
                phi_2ch = permute(phis{2}.big_mips(:, :, 3, channel_pair_counter, :, :, :, :), [1 5 6 7 8 2 3 4]);
                phi_2ch_weighted = phi_2ch;
                phi_2ch_weighted(:, :, :, :, :, 1) = phi_2ch_weighted(:, :, :, :, :, 1) .* single(permute(phis{2}.state_counters(:, channel_pair_counter, :, :, :, :), [1 3 4 5 6 2]));
                phi_2ch_weighted(:, :, :, :, :, 2) = phi_2ch_weighted(:, :, :, :, :, 2) .* single(permute(phis{2}.state_counters(:, channel_pair_counter, :, :, :, :), [1 3 4 5 6 2]));
                phi_2ch_summed = permute(sum(phi_2ch_weighted, 1), [2 3 4 5 6 1]); % Sum across states
                phi_2ch_mean = phi_2ch_summed / sum(phis{2}.state_counters(:, 1, 1, 1, 1, 4)); % Assumes all parameters have the same total number of states
                concept_phi.twoCh(match_counter, :, :, :, :, :) = phi_2ch_mean;
                
                % 4ch - get weighted average of small phi across states
                phi_4ch = permute(phis{4}.big_mips(:, :, concept, superset, :, :, :, :), [1 5 6 7 8 2 3 4]);
                phi_4ch_weighted = phi_4ch;
                phi_4ch_weighted(:, :, :, :, :, 1) = phi_4ch_weighted(:, :, :, :, :, 1) .* single(permute(phis{4}.state_counters(:, superset, :, :, :, :), [1 3 4 5 6 2]));
                phi_4ch_weighted(:, :, :, :, :, 2) = phi_4ch_weighted(:, :, :, :, :, 2) .* single(permute(phis{4}.state_counters(:, superset, :, :, :, :), [1 3 4 5 6 2]));
                phi_4ch_summed = permute(sum(phi_4ch_weighted, 1), [2 3 4 5 6 1]); % Sum across states
                phi_4ch_mean = phi_4ch_summed / sum(phis{4}.state_counters(:, 1, 1, 1, 1, 4)); % Assumes all parameters have the same total number of states
                concept_phi.fourCh(match_counter, :, :, :, :, :) = phi_4ch_mean;
                
                break;
            end
        end
        
        match_counter = match_counter + 1;
        
    end
end

%%

flies = (1:13);

figure;
%scatter(concept_phi.twoCh(:), concept_phi.fourCh(:));
x_vals = squeeze(mean(mean(concept_phi.twoCh(:, :, flies, 1, 4, 1), 2), 3));
y_vals = squeeze(mean(mean(concept_phi.fourCh(:, :, flies, 1, 4, 1), 2), 3));
scatter(x_vals, y_vals);

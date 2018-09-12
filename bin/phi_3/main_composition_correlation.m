%% Description

%{

For seeing if there's a relationship between two-channel phi and their
"corresponding" mechanism within 4-channel phi

%}

%% Setup

nChannels = [2 4];
tau = 1; % 1=4ms

%% Load relevant files

global_tpm = 0;

source_dir = 'results/';

phis = cell(4, 1);

% 2 channels
source_2ch = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels2_globalTPM' num2str(global_tpm) '.mat'];
tmp = load([source_dir source_2ch]);
phis{2} = tmp.phis{1};
phis{2}.concept_list_full = concept_list(2); % Relevant 2nd order concept is the 3rd one

% 4 channels
tic
disp('loading 4ch');
source_4ch = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM' num2str(global_tpm) '.mat'];
tmp = load([source_dir source_4ch]);
phis{4} = tmp.phis{1};
disp('loaded 4ch');
phis{4}.concept_list_full = concept_list(4);
toc

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
if length(tmp) < 8
    tmp(8) = 1;
end
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
        %channel_pair_location = channel_pair_location(2); Can compare first-order concepts as well
        % Get small-phi corresponding to the relevant concept
        for concept = 1 : length(phis{4}.concept_list_full)
            if isequal(phis{4}.concept_list_full{concept}, channel_pair_location)
                concept_2ch = 3; % The 2nd order concept in 2ch is always the third one (the last one)
                
                % 2ch - get weighted average of small phi across states
                phi_2ch = permute(phis{2}.big_mips(:, :, concept_2ch, channel_pair_counter, :, :, :, :), [1 5 6 7 8 2 3 4]);
                phi_2ch_weighted = phi_2ch;
                phi_2ch_weighted(:, :, :, :, :, 1) = phi_2ch_weighted(:, :, :, :, :, 1) .* single(permute(phis{2}.state_counters(:, channel_pair_counter, :, :, :, :), [1 3 4 5 6 2]));
                phi_2ch_weighted(:, :, :, :, :, 2) = phi_2ch_weighted(:, :, :, :, :, 2) .* single(permute(phis{2}.state_counters(:, channel_pair_counter, :, :, :, :), [1 3 4 5 6 2]));
                phi_2ch_summed = permute(sum(phi_2ch_weighted, 1), [2 3 4 5 6 1]); % Sum across states
                phi_2ch_mean = phi_2ch_summed / sum(phis{2}.state_counters(:, 1, 1, 1, 1, tau)); % Assumes all parameters have the same total number of states
                concept_phi.twoCh(match_counter, :, :, :, :, :) = phi_2ch_mean;
                
                % 4ch - get weighted average of small phi across states
                phi_4ch = permute(phis{4}.big_mips(:, :, concept, superset, :, :, :, :), [1 5 6 7 8 2 3 4]);
                phi_4ch_weighted = phi_4ch;
                phi_4ch_weighted(:, :, :, :, :, 1) = phi_4ch_weighted(:, :, :, :, :, 1) .* single(permute(phis{4}.state_counters(:, superset, :, :, :, :), [1 3 4 5 6 2]));
                phi_4ch_weighted(:, :, :, :, :, 2) = phi_4ch_weighted(:, :, :, :, :, 2) .* single(permute(phis{4}.state_counters(:, superset, :, :, :, :), [1 3 4 5 6 2]));
                phi_4ch_summed = permute(sum(phi_4ch_weighted, 1), [2 3 4 5 6 1]); % Sum across states
                phi_4ch_mean = phi_4ch_summed / sum(phis{4}.state_counters(:, 1, 1, 1, 1, tau)); % Assumes all parameters have the same total number of states
                concept_phi.fourCh(match_counter, :, :, :, :, :) = phi_4ch_mean;
                
                break;
            end
        end
        
        match_counter = match_counter + 1;
        
    end
end

%% Find portion of AB in ABCD for partitioned, unpartitioned
% Finds portion for each AB

% Find how many data entries we need
matches = 0;
for channel_pair_counter = 1 : size(phis{2}.channel_sets, 1)
    channel_pair = phis{2}.channel_sets(channel_pair_counter, :);
    
    supersets = find_supersets(channel_pair, phis{4}.channel_sets);
    matches = matches + length(supersets);
end

% Storage structure
tmp = size(phis{4}.big_mips);
if length(tmp) < 8
    tmp(8) = 1;
end
storage_size = [nchoosek(15, 2) tmp([5 6 7 8 2])];
portions = zeros(storage_size);

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
        %channel_pair_location = channel_pair_location(2); Can compare first-order concepts as well
        % Get small-phi corresponding to the relevant concept
        for concept = 1 : length(phis{4}.concept_list_full)
            if isequal(phis{4}.concept_list_full{concept}, channel_pair_location)
                concept_2ch = 3; % The 2nd order concept in 2ch is always the third one (the last one)
                
                % Find portion where concept was 0 in 4-channel set
                phi_4ch = permute(phis{4}.big_mips(:, :, concept, superset, :, :, :, :), [1 5 6 7 8 2 3 4]);
                phi_4ch = phi_4ch > 0;
                phi_4ch_weighted = zeros(size(phi_4ch));
                % Weight by number of occurrences of each state
                phi_4ch_weighted(:, :, :, :, :, 1) = phi_4ch(:, :, :, :, :, 1) .* single(permute(phis{4}.state_counters(:, superset, :, :, :, :), [1 3 4 5 6 2]));
                phi_4ch_weighted(:, :, :, :, :, 2) = phi_4ch(:, :, :, :, :, 2) .* single(permute(phis{4}.state_counters(:, superset, :, :, :, :), [1 3 4 5 6 2]));
                phi_4ch_summed = permute(sum(phi_4ch_weighted, 1), [2 3 4 5 6 1]); % Sum across states
                phi_4ch_mean = phi_4ch_summed / sum(phis{4}.state_counters(:, 1, 1, 1, 1, tau)); % Assumes all parameters have the same total number of states
                portions(channel_pair_counter, :, :, :, :, :) = portions(channel_pair_counter, :, :, :, :, :) + permute(phi_4ch_mean, [6 1 2 3 4 5]);
                
                break;
            end
        end
        
    end
    
    portions(channel_pair_counter, :, :, :, :, :) = portions(channel_pair_counter, :, :, :, :, :) ./ length(supersets);
    
end

figure;
titles = {'unpartitioned', 'partitioned'};
for partitioned = 1 : 2
    subplot(2, 1, partitioned);
    
    imagesc(squeeze(mean(mean(portions(:, :, :, :, 1, partitioned), 2), 3)));
    c = colorbar; title(c, 'portion AB present in ABCD');
    
    title(titles{partitioned});
    
    set(gca, 'XTick', [1 2]);
    set(gca, 'XTickLabel', {'wake', 'anest'});
end

ylabel('channel pair (AB)');
xlabel('condition');

%% Find portion of AB in ABCD for partitioned, unpartitioned
% Finds portion of 2nd order concepts which are 0, independent of reference 2nd order concept

portions = zeros(size(phis{4}.big_mips, 6), size(phis{4}.big_mips, 7), 2);

for fly = 1 : size(phis{4}.big_mips, 6)
    for condition = 1 : size(phis{4}.big_mips, 7)
        for partitioned = 1 : 2
            tmp = phis{4}.big_mips(:, partitioned, 5:5+6, :, :, fly, condition);
            portions(fly, condition, partitioned) = length(tmp(tmp == 0)) / numel(tmp);
        end
    end
end

figure;
titles = {'unpartitioned', 'partitioned'};
for partitioned = 1 : 2
    subplot(2, 1, partitioned);
    imagesc(portions(:, :, partitioned));
    colorbar;
    title(titles{partitioned});
end

%% Plot individually for each fly

flies = (1:13);

% figure;
% %scatter(concept_phi.twoCh(:), concept_phi.fourCh(:));
% x_vals = squeeze(mean(mean(concept_phi.twoCh(:, :, flies, 1, tau, 1), 2), 3));
% y_vals = squeeze(mean(mean(concept_phi.fourCh(:, :, flies, 1, tau, 1), 2), 3));
% scatter(x_vals, y_vals);

% Plot for each fly
figure;
condition = 1;
for fly = 1 : 13
    subplot(3, 6, fly);
    x_vals = squeeze(mean(mean(concept_phi.twoCh(:, :, fly, condition, tau, 1), 2), 3));
    y_vals = squeeze(mean(mean(concept_phi.fourCh(:, :, fly, condition, tau, 1), 2), 3));
    scatter(x_vals, y_vals);
    corr(x_vals, y_vals)
    
    axis equal;
    
    xlabel('2ch concept \phi'); ylabel('4ch concept \phi');
end

%% Plot values (mean across flies)
% Plot across flies
figure;
partition_titles = {'unpartitioned', 'partitioned'};
condition_colours = 'rb';
condition_shapes = 'ox';

% Averaged values across flies
flies = (1:13);
concept_phi_mean.twoCh = mean(mean(concept_phi.twoCh(:, :, flies, :, :, :), 2), 3);
concept_phi_mean.fourCh = mean(mean(concept_phi.fourCh(:, :, flies, :, :, :), 2), 3);
axis_lims = [min(concept_phi_mean.twoCh(:)) max(concept_phi_mean.twoCh(:)) min(concept_phi_mean.fourCh(:)) max(concept_phi_mean.fourCh(:))];

% Plot
for partitioned = 1 : 2
    subplot(1, size(concept_phi.twoCh, 4), partitioned);
    for condition = 1 : size(concept_phi.twoCh, 4)
        x_vals = squeeze(concept_phi_mean.twoCh(:, :, :, condition, tau, partitioned));
        y_vals = squeeze(concept_phi_mean.fourCh(:, :, :, condition, tau, partitioned));
        scatter(x_vals, y_vals, [condition_colours(condition) condition_shapes(condition)]); hold on;
        [r, p] = corr(x_vals, y_vals)
    end
    title(partition_titles{partitioned});
    
    axis(axis_lims);
    line(axis_lims(1:2), axis_lims(3:4), 'Color', 'k');
    
    xlabel('2ch concept \phi'); ylabel('4ch concept \phi');
end

legend('wake', 'anest', 'x=y', 'Location', 'southeast');

%% Plot difference (mean across flies)
% Plot across flies
figure;
partition_titles = {'unpartitioned', 'partitioned'};

% Averaged values across flies
flies = (1:13);
concept_phi_mean.twoCh = mean(mean(concept_phi.twoCh(:, :, flies, 1, :, :)-concept_phi.twoCh(:, :, flies, 2, :, :), 2), 3);
concept_phi_mean.fourCh = mean(mean(concept_phi.fourCh(:, :, flies, 1, :, :)-concept_phi.fourCh(:, :, flies, 2, :, :), 2), 3);
axis_lims = [min(concept_phi_mean.twoCh(:)) max(concept_phi_mean.twoCh(:)) min(concept_phi_mean.fourCh(:)) max(concept_phi_mean.fourCh(:))];

% Plot
for partitioned = 1 : 2
    subplot(1, size(concept_phi.twoCh, 4), partitioned);
    x_vals = squeeze(concept_phi_mean.twoCh(:, :, :, :, tau, partitioned));
    y_vals = squeeze(concept_phi_mean.fourCh(:, :, :, :, tau, partitioned));
    scatter(x_vals, y_vals, 'k.'); hold on;
    [r, p] = corr(x_vals, y_vals)
    title(partition_titles{partitioned});
    
    axis(axis_lims);
    line(axis_lims(1:2), axis_lims(3:4), 'Color', 'k');
    
    xlabel('2ch concept \phi'); ylabel('4ch concept \phi');
    
    ax = gca;
    ax.YAxis.Exponent = 0;
    ax.XAxis.Exponent = 0;
end

legend('wake-anest', 'x=y', 'Location', 'southeast');
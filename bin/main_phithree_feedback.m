%% DESCRIPTION

%{
This script tests the third hypothesis: phi-3 c->p cuts are more common under anaesthesia than under wakefulness

%}

%% SETUP
groupings.p = (1:6);
groupings.c = (9:14);

nChannels = 2; % This should always be 2, as we are looking at unidirectional c->p cuts

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    ];

%% LOAD
% disp('loading');
% 
% % Phi-3
% % .mips (holds the cuts) dimensions: states x sets x flies x conditions x taus
% load([data_directory data_filename '_phithree.mat']);
% 
% disp('loaded');

%% Reformat MIPs if necessary

% .mips should have dimensions corresponding to set, fly, condition, etc.
% If not, it will have size of length 2, with the first dimension being a combination of all parameters
if length(size(phis{1}.mips)) == 2
    
    % If partitions only consist of groups of one element (i.e. for 2 channels), then mips are stored in an array
    % We want to conver this to a cell
    if isa(phis{1}.mips, 'int64')
        phis{1}.mips = num2cell(phis{1}.mips);
    end
    
    % 4th dimension is nFlies (should be the same as specified during setup)
    [nStates, nSets, nFlies, nConditions, nTaus] = size(phis{1}.state_phis);
    
    % Looks like the partitions only splits into 2 groups
    % MIPs are calculated per state, not per trial
    % (only variance across trials is the state-weighting when averaging phi)
    mips_formatted = cell(nStates, nSets, nFlies, nConditions, nTaus);
    
    % Iterate through in same order as python3 script and reformat
    mip_counter = 1;
    for fly = 1 : nFlies
        for condition = 1 : nConditions
            for tau = 1 : nTaus
                for set = 1 : nSets
                    for state = 1 : nStates
                        cut{1} = phis{1}.mips(mip_counter, :);
                        mips_formatted(state, set, fly, condition, tau) = cut;
                        mip_counter = mip_counter + 1;
                    end
                end
            end
        end
    end
    
    % Replace the unformatted .mips with the new format
    phis{1}.mips = mips_formatted;
end

%% Relabel MIPs to have the actual channel label (instead of always 0, 1, 2, 3)

% for nChannels = 1 : length(phis)
%     for tau = 1 : size(phis{nChannels}.mips, 5)
%         for condition = 1 : size(phis{nChannels}.mips, 4)
%             for fly = 1 : size(phis{nChannels}.mips, 3)
%                 for set = 1 : size(phis{nChannels}.mips, 2)
%                     for state = 1 : size(phis{nChannels}.mips, 1)
%                         mip = phis{nChannels}.mips{state, set, fly, condition, tau};
%                         mip_relabelled = mip;
%                         for group = 1 : length(mip)
%                             for channel = 1 : length(mip{group})
%                                 mip_relabelled{group}(channel) = phis{nChannels}.channel_sets(set, mip{group}(channel)+1); % +1 because of python's 0 indexing
%                             end
%                         end
%                         phis{nChannels}.mips{state, set, fly, condition, tau} = mip_relabelled;
%                     end
%                 end
%             end
%         end
%     end
% end

%% Filter out sets which have a periphery and a centre channel

% Take only trials which contain a periphery channel and a central channel
[phis{1}.mips, phis{1}.state_counters, phis{1}.state_phis, phis{1}.channel_sets, phis{1}.phi_threes] = group_filter(phis{1}.mips, phis{1}.state_counters, phis{1}.state_phis, phis{1}.channel_sets, phis{1}.phi_threes, groupings);

for nChannels = 1 : length(phis)
    [phis{nChannels}.mips, phis{nChannels}.state_counters, phis{nChannels}.state_phis, phis{nChannels}.channel_sets, phis{nChannels}.phi_threes] = group_filter(phis{nChannels}.mips, phis{nChannels}.state_counters, phis{nChannels}.state_phis, phis{nChannels}.channel_sets, phis{nChannels}.phi_threes, groupings);
end

%% Portion of feedback cuts per trial, considering all MIPs

% Classify each cut in a trial as feedback or not
phis_unmoded.feedback = feedback_classify(phis{1}.mips, groupings);

% Count number of times states with a feedback cut occurred within the trial, for each channel set
phis_unmoded.feedback_counts = feedback_count(phis_unmoded.feedback, phis{1}.state_counters);

% Convert from counts to portions (percentage of cuts being feedback cuts)
phis_unmoded.feedback_portions = phis_unmoded.feedback_counts / sum(phis{1}.state_counters(:, 1, 1, 1, 1, 1));

% Average across sets and trials
%phis_unmoded.feedback_counts = squeeze(mean(mean(phis_unmoded.feedback_counts, 1), 2));

% Test for differing portions after averaging at each tau level
sigs = zeros(size(phis_unmoded.feedback_portions, 5), 1);
for tau = 1 : size(phis_unmoded.feedback_portions, 5)
    air = squeeze(mean(phis_unmoded.feedback_portions(:, :, :, 1, tau), 2));
    iso = squeeze(mean(phis_unmoded.feedback_portions(:, :, :, 2, tau), 2));
    
    [decision, sigs(tau)] = ttest(air(:), iso(:));
end

% Create barplot matrix (rows group bars together)
portion_plot = permute(squeeze(mean(mean(mean(phis_unmoded.feedback_portions, 1), 2), 3)), [2 1]);
portion_plot_err = permute(squeeze(std(mean(mean(phis_unmoded.feedback_portions, 1), 2), [], 3)), [2 1]) / sqrt(size(phis_unmoded.feedback_portions, 3));
figure;
bar((1:3)-0.15, portion_plot(:, 1), 0.25); hold on;
bar((1:3)+0.15, portion_plot(:, 2), 0.25, 'r');
errorbar((1:3)-0.15, portion_plot(:, 1), portion_plot_err(:, 1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
errorbar((1:3)+0.15, portion_plot(:, 2), portion_plot_err(:, 2), 'k', 'LineStyle', 'none', 'LineWidth', 1);
xticks((1:3));
xticklabels((phis{1}.taus));
xlabel('tau lag (ms)');
ylabel('c-/->p cut portion');


% Same plot format as raw phi values figure
labelled_subplot = 7;
plot_portions{1}.phis = squeeze(mean(mean(phis_unmoded.feedback_portions, 2), 3));
plot_portions{1}.phis_std = squeeze(std(mean(phis_unmoded.feedback_portions, 2), [], 3) / size(phis_unmoded.feedback_portions, 3));
plot_portions{1}.phis_delta = squeeze(mean(mean(phis_unmoded.feedback_portions(:, :, :, 1, :) - phis_unmoded.feedback_portions(:, :, :, 2, :), 2), 3));
plot_portions{1}.phis_delta_std = squeeze(std(mean(phis_unmoded.feedback_portions(:, :, :, 1, :) - phis_unmoded.feedback_portions(:, :, :, 2, :), 2), [], 3) / size(phis_unmoded.feedback_portions, 3));
plot_portions{1}.taus = phis{1}.taus;
plot_portions{1}.nChannels = phis{1}.nChannels;
plot_phis(plot_portions, [0 1], 1, 'portion', labelled_subplot);

%% Comparing phi of feedforward cuts to feedback cuts

% Filter out feedforward phis
ff_phis = phis{1}.state_phis(not(logical(phis_unmoded.feedback)));
% Filter out feedback phis
fb_phis = phis{1}.state_phis(logical(phis_unmoded.feedback));

% Filter out feedforward and feedback phis per fly (to find average and stderr)
ff_flies = zeros(size(phis{1}.state_phis, 3), size(phis{1}.state_phis, 4), size(phis{1}.state_phis, 5));
fb_flies = zeros(size(phis{1}.state_phis, 3), size(phis{1}.state_phis, 4), size(phis{1}.state_phis, 5));
for fly = 1 : size(phis{1}.state_phis, 3)
    for condition = 1 : size(phis{1}.state_phis, 4)
        for tau = 1 : size(phis{1}.state_phis, 5)
            state_phis = phis{1}.state_phis(:, :, fly, condition, tau);
            
            ff_filter = not(logical(phis_unmoded.feedback(:, :, fly, condition, tau)));
            fb_filter = logical(phis_unmoded.feedback(:, :, fly, condition, tau));
            
            ff_flies(fly, condition, tau) = mean(state_phis(ff_filter));
            fb_flies(fly, condition, tau) = mean(state_phis(fb_filter));
        end
    end
end

%% Portion of feedback cuts when considering only the most frequent MIP

% % Select most common phi-three MIP as the MIP for each trial
% disp('Choosing most frequent phi-3 MIPs');
% phis_moded = struct();
% [phis_moded.mips, phis_moded.mip_counts] = mip_mode(phis{1}.mips, phis{1}.state_counters);
% disp('Chosen');
% 
% % Classify that state's MIP as feedback or not
% phis_moded.feedback = feedback_classify(phis_moded.mips);
% 
% % Find portion of channel sets which have a c->p cut
% channel_sets = size(phis_moded.feedback, 1);
% phis_moded.feedback_counts = squeeze(sum(phis_moded.feedback, 1));

%% Likelihood of feedback cuts (during air) switching direction under anaesthesia
% Likelihood of feedforward cuts (during air) switching direction under anaesthesia
% Is this actually meaningful? Is 50/50 the default/chance probability?

% MIP - cut that makes the least difference to phi
% High feedback when awake = less likely that the MIP includes feedback cuts
% Normal feedback when anest = standard likelihood that the MIP includes feedback cuts
% Low frequency GC - we use tau=16 (longer lag) to compute TPM (what's the direct association?)
% Frequency range of reduced feedback = 1-10Hz, corresponding to down to 0.1s

%% Find only for sets which decreased significantly under iso
% % Find channels with positive delta phi
% 
% % Classify each cut in a trial as feedback or not
% phis_unmoded.feedback = feedback_classify(phis{1}.mips, groupings);
% 
% % Count number of times states with a feedback cut occurred within each trial, for each channel set
% phis_unmoded.feedback_counts = feedback_count(phis_unmoded.feedback, phis{1}.state_counters);
% 
% % Convert from counts to portions (percentage of cuts being feedback cuts)
% % Sum of state counters should be constant (2250)
% phis_unmoded.feedback_portions = phis_unmoded.feedback_counts / sum(phis{1}.state_counters(:, 1, 1, 1, 1, 1));
% 
% % Filter for channel sets and trials where air > iso
% % Because different flies will likely have a different number of sets, the result should be a cell array, with a cell per fly
% phis_unmoded.feedback_portions_filtered = phi_delta_filter(phis_unmoded.feedback_portions, phis{1}.phi_threes);
% 
% % Get average portions per fly and tau
% % i.e. average across trials and filtered sets
% phis_unmoded.feedback_portions_filtered_trialmeans = zeros(length(phis_unmoded.feedback_portions_filtered), size(phis_unmoded.feedback_portions, 4), size(phis_unmoded.feedback_portions, 5));
% for fly = 1 : length(phis_unmoded.feedback_portions_filtered)
%     for condition = 1 : size(phis_unmoded.feedback_portions, 4)
%         for tau = 1 : size(phis_unmoded.feedback_portions, 5)
%             phis_unmoded.feedback_portions_filtered_trialmeans(fly, condition, tau) = mean(mean(phis_unmoded.feedback_portions_filtered{fly, tau}(:, :, condition), 2), 1);
%         end
%     end
% end
% 
% % Test for differing portions at each tau level
% sigs = zeros(size(phis_unmoded.feedback_portions_filtered_trialmeans, 3), 1);
% for tau = 1 : size(phis_unmoded.feedback_portions, 5)
%     air = phis_unmoded.feedback_portions_filtered_trialmeans(:, 1, tau);
%     iso = phis_unmoded.feedback_portions_filtered_trialmeans(:, 2, tau);
%     
%     [decision, sigs(tau)] = ttest(air(:), iso(:));
% end
% 
% % Create barplot matrix (rows group bars together)
% portion_plot = permute(squeeze(mean(phis_unmoded.feedback_portions_filtered_trialmeans, 1)), [2 1]);
% portion_plot_err = permute(squeeze(std(phis_unmoded.feedback_portions_filtered_trialmeans, [], 1)), [2 1]) / sqrt(size(phis_unmoded.feedback_portions_filtered_trialmeans, 1));
% figure;
% bar((1:3)-0.15, portion_plot(:, 1), 0.25); hold on;
% bar((1:3)+0.15, portion_plot(:, 2), 0.25, 'r');
% errorbar((1:3)-0.15, portion_plot(:, 1), portion_plot_err(:, 1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
% errorbar((1:3)+0.15, portion_plot(:, 2), portion_plot_err(:, 2), 'k', 'LineStyle', 'none', 'LineWidth', 1);
% xticks((1:3));
% xticklabels((phis{1}.taus));
% xlabel('tau lag (ms)');
% ylabel('c-/->p cut portion');
% 
% function [filtered] = phi_delta_filter(portions, phis)
% % Select only trials and sets which have higher phi under air
% % Averages differences across trials, sets which on average have air>iso are kept
% % Because different flies and taus will likely have different filtered sets and trials, a cell is returned
% % where each cell corresponds to a fly and tau
% 
% filtered = cell(size(phis, 3), size(phis, 5));
% for fly = 1 : size(phis, 3)
%     for tau = 1 : size(phis, 5)
%         fly_phis = squeeze(mean(phis, 2));
%         valid_delta = fly_phis(:, fly, 1, tau) > fly_phis(:, fly, 2, tau);
%         filtered{fly, tau} = squeeze(portions(valid_delta, :, fly, :, tau));
%     end
% end
% 
% end


%% Function: grouping filter

function [filtered_mips, filtered_counters, filtered_state_phis, filtered_sets, filtered_phis] = group_filter(mips, state_counters, state_phis, channel_sets, phis, groupings)
% Filters for mips which contain channels as specified by groupings
%

%%%%%%%%%%%%
% First channel is in periphery and second is in centre, or vice versa
% However, if channel_sets is sorted, the second case never occurs (so it shouldn't be sorted)
% This is the 2 channel condition
filter_indices =...
    ((channel_sets(:, 1) >= groupings.p(1) & channel_sets(:, 1) <= groupings.p(end)) &...
    (channel_sets(:, 2) >= groupings.c(1) & channel_sets(:, 2) <= groupings.c(end))) |...
    ((channel_sets(:, 1) >= groupings.c(1) & channel_sets(:, 1) <= groupings.c(end)) &...
    (channel_sets(:, 2) >= groupings.p(1) & channel_sets(:, 2) <= groupings.p(end)));
%%%%%%%%%%%

%%%%%%%%%%%
% % This is for selecting sets which are built of the grouping channels (with no other restrictions)
% filter_indices = zeros(size(channel_sets, 1), 1);
% % For each channel set, check if set is contains only the channels in grouping
% for set = 1 : size(channel_sets, 1)
%     within_grouping = true;
%     for channel = 1 : size(channel_sets, 2)
%         % Check if channel is peripheral OR central
%         % If it is not, then within_grouping should become 0
%         % If it is, and within_grouping is 0, then it should stay 0
%         within_grouping = min([within_grouping (any(groupings.p==channel_sets(set, channel)) || any(groupings.c==channel_sets(set, channel)))]);
%     end
%     filter_indices(set) = within_grouping;
% end
%%%%%%%%%%%%%%

%%%%%%%%%%%%%
% This is for selecting sets where each partition group consists of either ONLY peripheral channels, or ONLY central channels
% Actually, not doable, as this depends on the MIP, not the channel set (and the MIP varies per state within a single channel set)
%%%%%%%%%%%%%

filter_indices = logical(filter_indices);

filtered_sets = channel_sets(filter_indices, :);

filtered_mips = mips(:, filter_indices, :, :, :);

filtered_state_phis = state_phis(:, filter_indices, :, :, :);

filtered_counters = state_counters(:, filter_indices, :, :, :, :);

filtered_phis = phis(filter_indices, :, :, :, :);

end

%% Function: count number of feedback cuts

function [counts] = feedback_count(classified, frequencies)
% This is for multiple mips per trial (not selecting the most frequent mip as the mip)
% Inputs:
%   classified = logical array, corresponding to feedback-classified mips
% Outputs:
%   portions = 

portion_dimensions = size(frequencies);

counts = zeros(portion_dimensions(2:end));
for tau = 1 : size(frequencies, 6)
    for condition = 1 : size(frequencies, 5)
        for fly = 1 : size(frequencies, 4)
            for trial = 1 : size(frequencies, 3)
                for set = 1 : size(frequencies, 2)
                    counts(set, trial, fly, condition, tau) = sum(squeeze(frequencies(logical(classified(:, set, fly, condition, tau)), set, trial, fly, condition, tau)));
                end
            end
        end
    end
end

end

%% Function: determine if cut is feedback or feedfoward

function [feedback] = feedback_cut(partition_cut, groupings)
% Returns 1 if cut is c -> p (i.e. first element is larger than second; assumes larger elements are more central)
%
% A cut is defined as feedback if it is made of more feedback cuts than feedforward cuts
% e.g. for A B C D (A is most peripheral), if the cut is AB-/->CD, then it is made of
% the cuts A-/->C, A-/->D, B-/->C, and B-/->D, so 4 feedback cuts and no feedforward
% e.g. for A B C D (A is most peripheral), if the cut is AD-/->BC, then it is made of
% the cuts A-/->B, A-/->C, D-/->B, and D-/->C, so 2 feedback cuts, and 2 feedforward
% Due to channel labelling convention, D>C>B>A
%
% Inputs:
%   cut = cell array with length 2, each cell holds a single numeric
%   cut = cell array representing a partition cut (i.e. with length 2); the cut is from
%       the channels in the first cell to the channels in the second cell
% Outputs:
%   feedback = 1 if the cut is to a feedback connection, 0 otherwise (if equal feedback and feedforward, 0)
cuts_fb = 0;
cuts_ff = 0;

for from = 1 : length(partition_cut{1}) % for each cut source
    for to = 1 : length(partition_cut{2}) % for each cut end
        if any(groupings.c==partition_cut{1}(from)) && any(groupings.p==partition_cut{2}(to)) % feedback if source is central and target is peripheral
            cuts_fb = cuts_fb + 1;
        else
            cuts_ff = cuts_ff + 1;
        end
    end
end

if cuts_fb > cuts_ff
    feedback = 1;
else
    feedback = 0;
end

end

%% Function: classify cuts as feedback or not

function [classified] = feedback_classify(mips, groupings)
% Classifies each mip/cut in mips as feedback (1) or not (0)
%
% Inputs:
%   mips = cell array holding cuts
% Outputs:
%   classified = numeric array (logical) with 1s where a cut corresponds to feedback

classified = zeros(size(mips));

for cut = 1 : numel(mips)
    classified(cut) = feedback_cut(mips{cut}, groupings);
end

end

%% Function: concatenate all portions and plot on a single axis

function [] = plot_phis(phis, y_limits, errorbars, ylabel_text, labelled_subplot)
% Plots all phis (for all nChannels used) on a single axis
%
% Inputs:
%   phis: cell array, each cell holds .phis, a matrix with dimensions
%       (channel sets x conditions x taus); Average across trials and average across
%       flies or or filter for a fly before inputting to this function, .phis_std, .phis_delta, and .phis_delta_std
% Outputs: None

plot_titles{1} = 'air';
plot_titles{2} = 'iso';
plot_titles{3} = 'delta';
tau_labels = phis{1}.taus;

phis_all = [];
phis_all_stds = [];
phis_all_delta = [];
phis_all_delta_std = [];
channel_ticks = []; xtick_counter = 1;
channel_labels = [];
% Concatenate results
for nChannels_counter = 1 : length(phis)
    phis_all = cat(1, phis_all, phis{nChannels_counter}.phis);
    phis_all_stds = cat(1, phis_all_stds, phis{nChannels_counter}.phis_std);
    phis_all_delta = cat(1, phis_all_delta, phis{nChannels_counter}.phis_delta);
    phis_all_delta_std = cat(1, phis_all_delta_std, phis{nChannels_counter}.phis_delta_std);
    channel_ticks = cat(1, channel_ticks, xtick_counter);
    channel_labels = cat(1, channel_labels, phis{nChannels_counter}.nChannels);
    xtick_counter = xtick_counter + size(phis{nChannels_counter}.phis, 1);
end

% Plot for each tau, for each condition, averaging across trials and flies
figure;
subplot_counter = 1;
for tau = 1 : size(phis_all, 3)
    for condition = 1 : size(phis_all, 2)
        subplot(size(phis_all, 3), size(phis_all, 2) + 1, subplot_counter);
        
        %barwitherr(squeeze(std(mean(phis_all(:, :, :, condition, tau), 2), [], 3)), squeeze(mean(mean(phis_all(:, :, :, condition, tau), 2), 3)), 'LineWidth', 0.001);
        
        bar(phis_all(:, condition, tau));
        if subplot_counter <= size(phis_all, 2) + 1
            title(plot_titles{subplot_counter});
        end
        if errorbars == 1
            hold on;
            errorbar((1:size(phis_all, 1)), phis_all(:, condition, tau), zeros(size(phis_all, 1), 1), phis_all_stds(:, condition, tau), 'k', 'LineStyle', 'none', 'LineWidth', 0.1, 'CapSize', 0);
        end
        axis([-1 size(phis_all, 1)+2 y_limits]);
        xticks(channel_ticks); xticklabels(channel_labels);
        
        if subplot_counter == labelled_subplot
            xlabel('channel set');
            ylabel([ylabel_text ' tau=' num2str(tau_labels(tau))]);
        elseif mod(subplot_counter, size(phis_all, 3)) == 1
            ylabel(['tau=' num2str(tau_labels(tau))]);
        end
        
        subplot_counter = subplot_counter + 1;
    end
    
    % Plot delta
    subplot(size(phis_all, 3), size(phis_all, 2) + 1, subplot_counter);
    bar(phis_all_delta(:, tau));
    if subplot_counter <= size(phis_all, 2) + 1
        title(plot_titles{subplot_counter});
    end
    if errorbars == 1
        hold on;
        errorbar((1:size(phis_all, 1)), phis_all_delta(:, tau), zeros(size(phis_all, 1), 1), phis_all_delta_std(:, tau), 'k', 'LineStyle', 'none', 'LineWidth', 0.1, 'CapSize', 0);
    end
    axis([-1 size(phis_all, 1)+2 -0.4 0.4]);
    xticks(channel_ticks); xticklabels(channel_labels);
    
    if subplot_counter + 1 - size(phis_all, 3) == labelled_subplot
        ylabel(['delta ' ylabel_text]);
    end
    
    subplot_counter = subplot_counter + 1;
end

end
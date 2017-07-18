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
% If not, it will have size 2, with the first dimension being a combination of all parameters
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

%% Filter out sets which have a periphery and a centre channel

% Take only trials which contain a periphery channel and a central channel
[phis{1}.mips, phis{1}.state_counters, phis{1}.channel_sets, phis{1}.phi_threes] = group_filter(phis{1}.mips, phis{1}.state_counters, phis{1}.channel_sets, phis{1}.phi_threes, groupings);

%% Portion of feedback cuts per trial, considering all MIPs

% % Classify each cut in a trial as feedback or not
% phis_unmoded.feedback = feedback_classify(phis{1}.mips);
% 
% % Count number of times states with a feedback cut occurred within the trial, for each channel set
% phis_unmoded.feedback_counts = feedback_count(phis_unmoded.feedback, phis{1}.state_counters);
% 
% % Convert from counts to portions (percentage of cuts being feedback cuts)
% phis_unmoded.feedback_portions = phis_unmoded.feedback_counts / sum(phis{1}.state_counters(:, 1, 1, 1, 1, 1));
% 
% % Average across sets and trials
% %phis_unmoded.feedback_counts = squeeze(mean(mean(phis_unmoded.feedback_counts, 1), 2));
% 
% % Test for differing portions after averaging at each tau level
% sigs = zeros(size(phis_unmoded.feedback_portions, 5), 1);
% for tau = 1 : size(phis_unmoded.feedback_portions, 5)
%     air = squeeze(mean(phis_unmoded.feedback_portions(:, :, :, 1, tau), 2));
%     iso = squeeze(mean(phis_unmoded.feedback_portions(:, :, :, 2, tau), 2));
%     
%     [decision, sigs(tau)] = ttest(air(:), iso(:));
% end
% 
% % Create barplot matrix (rows group bars together)
% portion_plot = permute(squeeze(mean(mean(mean(phis_unmoded.feedback_portions, 1), 2), 3)), [2 1]);
% portion_plot_err = permute(squeeze(std(mean(mean(phis_unmoded.feedback_portions, 1), 2), [], 3)), [2 1]) / sqrt(size(phis_unmoded.feedback_portions, 3));
% figure;
% bar((1:3)-0.15, portion_plot(:, 1), 0.25); hold on;
% bar((1:3)+0.15, portion_plot(:, 2), 0.25, 'r');
% errorbar((1:3)-0.15, portion_plot(:, 1), portion_plot_err(:, 1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
% errorbar((1:3)+0.15, portion_plot(:, 2), portion_plot_err(:, 2), 'k', 'LineStyle', 'none', 'LineWidth', 1);
% xticks((1:3));
% xticklabels((phis{1}.taus));
% xlabel('tau lag (ms)');
% ylabel('c-/->p cut portion');

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
% Find channels with positive delta phi

% Classify each cut in a trial as feedback or not
phis_unmoded.feedback = feedback_classify(phis{1}.mips);

% Count number of times states with a feedback cut occurred within each trial, for each channel set
phis_unmoded.feedback_counts = feedback_count(phis_unmoded.feedback, phis{1}.state_counters);

% Convert from counts to portions (percentage of cuts being feedback cuts)
% Sum of state counters should be constant (2250)
phis_unmoded.feedback_portions = phis_unmoded.feedback_counts / sum(phis{1}.state_counters(:, 1, 1, 1, 1, 1));

% Filter for channel sets and trials where air > iso
% Because different flies will likely have a different number of sets, the result should be a cell array, with a cell per fly
phis_unmoded.feedback_portions_filtered = phi_delta_filter(phis_unmoded.feedback_portions, phis{1}.phi_threes);

% Get average portions per fly and tau
% i.e. average across trials and filtered sets
phis_unmoded.feedback_portions_filtered_trialmeans = zeros(length(phis_unmoded.feedback_portions_filtered), size(phis_unmoded.feedback_portions, 4), size(phis_unmoded.feedback_portions, 5));
for fly = 1 : length(phis_unmoded.feedback_portions_filtered)
    for condition = 1 : size(phis_unmoded.feedback_portions, 4)
        for tau = 1 : size(phis_unmoded.feedback_portions, 5)
            phis_unmoded.feedback_portions_filtered_trialmeans(fly, condition, tau) = mean(mean(phis_unmoded.feedback_portions_filtered{fly, tau}(:, :, condition), 2), 1);
        end
    end
end

% Test for differing portions at each tau level
sigs = zeros(size(phis_unmoded.feedback_portions_filtered_trialmeans, 3), 1);
for tau = 1 : size(phis_unmoded.feedback_portions, 5)
    air = phis_unmoded.feedback_portions_filtered_trialmeans(:, 1, tau);
    iso = phis_unmoded.feedback_portions_filtered_trialmeans(:, 2, tau);
    
    [decision, sigs(tau)] = ttest(air(:), iso(:));
end

% Create barplot matrix (rows group bars together)
portion_plot = permute(squeeze(mean(phis_unmoded.feedback_portions_filtered_trialmeans, 1)), [2 1]);
portion_plot_err = permute(squeeze(std(phis_unmoded.feedback_portions_filtered_trialmeans, [], 1)), [2 1]) / sqrt(size(phis_unmoded.feedback_portions_filtered_trialmeans, 1));
figure;
bar((1:3)-0.15, portion_plot(:, 1), 0.25); hold on;
bar((1:3)+0.15, portion_plot(:, 2), 0.25, 'r');
errorbar((1:3)-0.15, portion_plot(:, 1), portion_plot_err(:, 1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
errorbar((1:3)+0.15, portion_plot(:, 2), portion_plot_err(:, 2), 'k', 'LineStyle', 'none', 'LineWidth', 1);
xticks((1:3));
xticklabels((phis{1}.taus));
xlabel('tau lag (ms)');
ylabel('c-/->p cut portion');

function [filtered] = phi_delta_filter(portions, phis)
% Select only trials and sets which have higher phi under air
% Averages differences across trials, sets which on average have air>iso are kept
% Because different flies and taus will likely have different filtered sets and trials, a cell is returned
% where each cell corresponds to a fly and tau

filtered = cell(size(phis, 3), size(phis, 5));
for fly = 1 : size(phis, 3)
    for tau = 1 : size(phis, 5)
        fly_phis = squeeze(mean(phis, 2));
        valid_delta = fly_phis(:, fly, 1, tau) > fly_phis(:, fly, 2, tau);
        filtered{fly, tau} = squeeze(portions(valid_delta, :, fly, :, tau));
    end
end

end


%% Function: grouping filter

function [filtered_mips, filtered_counters, filtered_sets, filtered_phis] = group_filter(mips, state_counters, channel_sets, phis, groupings)
% Filters for mips which contain channels as specified by groupings
%

% First channel is in periphery and second is in centre, or vice versa
% However, if channel_sets is sorted, the second case never occurs
filter_indices =...
    ((channel_sets(:, 1) >= groupings.p(1) & channel_sets(:, 1) <= groupings.p(end)) &...
    (channel_sets(:, 2) >= groupings.c(1) & channel_sets(:, 2) <= groupings.c(end))) |...
    ((channel_sets(:, 1) >= groupings.c(1) & channel_sets(:, 1) <= groupings.c(end)) &...
    (channel_sets(:, 2) >= groupings.p(1) & channel_sets(:, 2) <= groupings.p(end)));

filtered_sets = channel_sets(filter_indices, :);

filtered_mips = mips(:, filter_indices, :, :, :);

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

function [feedback] = feedback_cut(cut)
% Returns 1 if cut is c -> p (i.e. first element is larger than second; assumes larger elements are more central)
%
% Inputs:
%   cut = cell array with length 2, each cell holds a single numeric
% Outputs:
%   feedback = 1 if the cut is to a feedback connection, 0 otherwise

feedback = 0;

if cut{1} > cut{2}
    feedback = 1;
end

end

%% Function: classify cuts as feedback or not

function [classified] = feedback_classify(mips)
% Classifies each mip/cut in mips as feedback (1) or not (0)
%
% Inputs:
%   mips = cell array holding cuts
% Outputs:
%   classified = numeric array (logical) with 1s where a cut corresponds to feedback

classified = zeros(size(mips));

for cut = 1 : numel(mips)
    classified(cut) = feedback_cut(mips{cut});
end

end
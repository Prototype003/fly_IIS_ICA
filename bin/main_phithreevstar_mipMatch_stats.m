%% DESCRIPTION

%{
This script tests the second hypothesis: phi-3 MIP cuts are equivalent to phi-star MIPs

Note that phi-3 only bipartitions

Chance level of match for 3 channels:
A BC
AB C
AC B
3/9 = 1/3 (for both bipartition filtering and no filtering)
Bell(3) = 5

Chance level of match for 4 channels:

A BCD
AB CD
AC BD
AD BC
ABC D
ABD C
ACD B
7/42    (not 1/7 if not filtering for bipartitions)
1/7 if filtering for bipartitions only
Bell(4) = 15


MIP* == MIP3, per set, trial, ...

A:
MIP* per set, trial, ...
2250 MIP3 per set, trial, ...
%match per set, trial, ...

B:
MIP* per set, trial, ...
MIP3 per set, trial, ...
1/0 per set, trial, ...
%}

%% SETUP

chance_levels = [1 1/3 1/7];% [1 1/3 7/42]; (7/42 is wrong....)

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_phistar.mat'
    ];

results_directory = 'analysis_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_mipComparison.mat'
    ];

%% LOAD
disp('loading');

% Match results
load([results_directory results_filename]);

% Phi-star results (for exclusion of certain MIPs)
load([data_directory data_filename]);

matches_per_trial_bkp = matches_per_trial;
matches_bkp = matches;

disp('loaded')

%% Filter results

% Exclude results from MIPs which were not bipartitions (because these never match with phi-3 bipartitioning)
% In order to exclude MIPs which are not bipartitions, we sum the number of occurrences together, ignoring nonbipartitions
% Then average afterwards, correcting for the total number of samples
% this works because the number of samples is constant, so avg(total across trials) = avg(avg per trial)

matches_old = matches;
matches = cell(length(matches_old), 1); % update matches to have a matching format with matches_per_trial
for nChannels_counter = 1 : length(matches_per_trial)
    [matches_per_trial{nChannels_counter}.matches_sum, matches_per_trial{nChannels_counter}.trial_totals] = sum_bipartitions_only(phis{nChannels_counter}.mips, matches_per_trial{nChannels_counter}.matches, matches_per_trial{nChannels_counter}.trial_total);
    matches{nChannels_counter}.matches = matches_old{nChannels_counter};
    [matches{nChannels_counter}.matches_sum, matches{nChannels_counter}.trial_totals] = sum_bipartitions_only(phis{nChannels_counter}.mips, matches{nChannels_counter}.matches, 1);
end

%% Overall matching percentage (without filtering for bipartition MIPs)

% match_portions_per_trial = cell(length(matches_per_trial), 1);
% for nChannels_counter = 1 : length(matches_per_trial)
%     
%     % Get portion of matches per trial
%     match_portions_per_trial{nChannels_counter} = matches_per_trial{nChannels_counter}.matches / matches_per_trial{nChannels_counter}.trial_total;
%     
%     % Average across channel sets and trials and flies
%     match_portions_per_trial{nChannels_counter} = mean(mean(mean(match_portions_per_trial{nChannels_counter}, 1), 2), 3);
%     match_portions_per_trial{nChannels_counter} = squeeze(match_portions_per_trial{nChannels_counter});
%     
%     
% end

%% Overall matching percentage (after filtering for bipartition MIPs)

match_portions_per_trial_filtered = cell(length(matches_per_trial), 1);
match_portions_filtered = cell(length(matches), 1);
for nChannels_counter = 1 : length(matches_per_trial)
    
    % Get portions of matches per trial
    match_portions_per_trial_filtered{nChannels_counter} = matches_per_trial{nChannels_counter}.matches_sum ./ matches_per_trial{nChannels_counter}.trial_totals;
    
    % Average across trials
    match_portions_per_trial_filtered{nChannels_counter} = squeeze(mean(match_portions_per_trial_filtered{nChannels_counter}, 1));
    
    
    match_portions_filtered{nChannels_counter} = matches{nChannels_counter}.matches_sum ./ matches{nChannels_counter}.trial_totals;
    match_portions_filtered{nChannels_counter} = squeeze(mean(match_portions_filtered{nChannels_counter}, 1));
end

% Plot
figure;
for nChannels_counter = 2 : length(matches_per_trial)
    subplot(1, length(matches_per_trial)-1, nChannels_counter-1);
    match_portion = permute(squeeze(mean(match_portions_per_trial_filtered{nChannels_counter}, 1)), [2 1]);
    match_portion_err = permute(squeeze(std(match_portions_per_trial_filtered{nChannels_counter}, [], 1)), [2, 1]) / sqrt(size(match_portions_per_trial_filtered{nChannels_counter}, 1));
    
    bar((1:3)-0.15, match_portion(:, 1), 0.25); hold on;
    bar((1:3)+0.15, match_portion(:, 2), 0.25, 'r');
    errorbar((1:3)-0.15, match_portion(:, 1), match_portion_err(:, 1), 'k', 'LineStyle', 'none', 'LineWidth', 1);
    errorbar((1:3)+0.15, match_portion(:, 2), match_portion_err(:, 2), 'k', 'LineStyle', 'none', 'LineWidth', 1);
    line([0 length(matches_per_trial)+1], [chance_levels(nChannels_counter) chance_levels(nChannels_counter)], 'Color', 'black', 'LineStyle', ':');
    
    title([num2str(phis{nChannels_counter}.nChannels) ' channels']);
    xticks((1:3));
    xticklabels((phis{1}.taus));
    axis([0 4 0 0.5]);
    
    if nChannels_counter == 2
        xlabel('tau lag (ms)');
        ylabel('match portion');
    end
end

%% Function: Remove results which do not correspond to a phi-3 bipartition
% phi-star partitions which are not bipartitions will never match phi-3 bipartitions

function [bipartition_sum, total_correction] = sum_bipartitions_only(mips, matches, samples_per_trial)
% Adds channel sets which correspond to bipartitions
% Sums across sets
%
% Inputs:
%
% Outputs:
%

sum_dimensions = size(mips);

bipartition_sum = zeros(sum_dimensions(2:end));
total_correction = zeros(sum_dimensions(2:end));

for tau = 1 : size(mips, 5)
    for condition = 1 : size(mips, 4)
        for fly = 1 : size(mips, 3)
            for trial = 1 : size(mips, 2)
                for set = size(mips, 1) : -1 : 1
                    if length(mips{set, trial, fly, condition, tau}) == 2
                        bipartition_sum(trial, fly, condition, tau) = bipartition_sum(trial, fly, condition, tau) + matches(set, trial, fly, condition, tau);
                        total_correction(trial, fly, condition, tau) = total_correction(trial, fly, condition, tau) + samples_per_trial;
                    end
                end
            end
        end
    end
end

end
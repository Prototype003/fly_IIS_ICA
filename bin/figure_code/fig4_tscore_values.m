%% Description

%{

Obtains values for plotting via fig4_tscore.m

%}

%% Setup

measure = 'phi_three';
tau = 1; % 1 = 4ms; 2 = 8ms; 3 = 16ms
if tau == 1
    tau_string = '4';
elseif tau == 2
    tau_string = '8';
elseif tau == 3
    tau_string = '16';
end

freq_range_w = (1:42); %(1:83); % corresponding to ~5Hz and ~10Hz, check the 'frequencies' vector
freq_range_a = (1:329); %(1:329)=0-5Hz; There are more frequency bins for the single large trial
freq_range_string = '0-5Hz'; %'0-10Hz';

fontsize = 11; % Used for drawing label letters

bin_location = '../';
addpath(bin_location);

results_directory = [bin_location 'workspace_results/'];

out_file = 'fig4_rankscore.mat';

%% Data for boxplots (within)

% measure_accuracies - this held classification accuracy, but now will hold t-scores

% Power
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification.mat';
load([results_directory results_filename]);
powers = permute(mean(powers(freq_range_w, :, :, :, :), 1), [2 3 5 4 1]); % mean across frequency range; trials x channels x flies x conditions
scores = zeros(size(powers, 3), size(powers, 2)); % flies x channels
for channel = 1 : size(powers, 2)
    for fly = 1 : size(powers, 3)
        % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
        [h, p, ci, stats] = ttest2(powers(:, channel, fly, 1), powers(:, channel, fly, 2));
        [p, h, stats] = ranksum(powers(:, channel, fly, 1), powers(:, channel, fly, 2));
        scores(fly, channel) = stats.ranksum; % stats.ranksum; stats.tstat;
    end
end
% Average scores across channels
scores_mean = mean(scores, 2);
% Add to boxplot data structure
measure_accuracies_w = scores_mean;
measure_groups_w = zeros(size(measure_accuracies_w)) + 1;



% Coherence
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';
load([results_directory results_filename]);
coherencies = permute(mean(coherencies(freq_range_w, :, :, :, :), 1), [2 3 5 4 1]); % mean across frequency range; trials x ch-pairs x flies x conditions
scores = zeros(size(coherencies, 3), size(coherencies, 2)); % flies x ch-pairs
for pair = 1 : size(coherencies, 2)
    for fly = 1 : size(coherencies, 3)
        % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
        [h, p, ci, stats] = ttest2(coherencies(:, pair, fly, 1), coherencies(:, pair, fly, 2));
        [p, h, stats] = ranksum(coherencies(:, pair, fly, 1), coherencies(:, pair, fly, 2));
        scores(fly, pair) = stats.ranksum; % stats.ranksum; stats.tstat;
    end
end
% Average scores across pairs
scores_mean = mean(scores, 2);
% Add to boxplot data structure
measure_accuracies_w = [measure_accuracies_w; scores_mean];
measure_groups_w = [measure_groups_w; zeros(size(scores_mean)) + 2];



% Phi-three
[phis, measure_string] = phi_load('phi_three', 0, bin_location);
group_counter = 3;
for nChannels_counter = 1 : length(phis)
    values = phis{nChannels_counter}.phis(:, :, :, :, tau);
    phis{nChannels_counter}.tscores = zeros(size(values, 3), size(values, 1)); % flies x ch-sets
    for set_counter = 1 : size(values, 1)
        for fly = 1 : size(values, 3)
            % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
            [h, p, ci, stats] = ttest2(values(set_counter, :, fly, 1), values(set_counter, :, fly, 2));
            [p, h, stats] = ranksum(values(set_counter, :, fly, 1), values(set_counter, :, fly, 2));
            phis{nChannels_counter}.tscores(fly, set_counter) = stats.ranksum; %stats.ranksum; stats.tstat;
        end
    end
    % Average scores across sets
    scores_mean = mean(phis{nChannels_counter}.tscores, 2);
    % Add to boxplot data structure
    measure_accuracies_w = [measure_accuracies_w; scores_mean];
    measure_groups_w = [measure_groups_w; zeros(size(scores_mean)) + group_counter];
    group_counter = group_counter + 1;
end
% Store for other plots (store in accuracies - in order to reuse code)
if strcmp(measure, 'phi_three')
    accuracies_w = cell(size(phis));
    for nChannels_counter = 1 : length(phis)
        accuracies_w{nChannels_counter}.accuracies = permute(phis{nChannels_counter}.tscores, [2 1]);
        accuracies_w{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
    end
end



% Phi-star
[phis, measure_string] = phi_load('phi_star', 0, bin_location);
for nChannels_counter = 1 : length(phis)
    values = phis{nChannels_counter}.phis(:, :, :, :, tau);
    phis{nChannels_counter}.tscores = zeros(size(values, 3), size(values, 1)); % flies x ch-sets
    for set_counter = 1 : size(values, 1)
        for fly = 1 : size(values, 3)
            % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
            [h, p, ci, stats] = ttest2(values(set_counter, :, fly, 1), values(set_counter, :, fly, 2));
            [p, h, stats] = ranksum(values(set_counter, :, fly, 1), values(set_counter, :, fly, 2));
            phis{nChannels_counter}.tscores(fly, set_counter) = stats.ranksum; %stats.ranksum; stats.tstat;
        end
    end
    % Average scores across sets
    scores_mean = mean(phis{nChannels_counter}.tscores, 2);
    % Add to boxplot data structure
    measure_accuracies_w = [measure_accuracies_w; scores_mean];
    measure_groups_w = [measure_groups_w; zeros(size(scores_mean)) + group_counter];
    group_counter = group_counter + 1;
end
% Store for other plots (store in accuracies - in order to reuse code)
if strcmp(measure, 'phi_star')
    accuracies_w = cell(size(phis));
    for nChannels_counter = 1 : length(phis)
        accuracies_w{nChannels_counter}.accuracies = permute(phis{nChannels_counter}.tscores, [2 1]);
        accuracies_w{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
    end
end

%% Data for boxplots (across)

% Power
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across1.mat';
load([results_directory results_filename]);
powers = permute(mean(powers(freq_range_a, :, :, :, :), 1), [5 3 4 1 2]); % mean across frequency range; flies x channels x conditions
scores = zeros(size(powers, 2), 1); % channels
for channel = 1 : size(powers, 2)
    % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
    [h, p, ci, stats] = ttest(powers(:, channel, 1), powers(:, channel, 2));
    [p, h, stats] = signrank(powers(:, channel, 1), powers(:, channel, 2));
    scores(channel) = stats.signedrank; % stats.signedrank; stats.tstat;
end
% Add to boxplot data structure
measure_accuracies_a = scores;
measure_groups_a = zeros(size(measure_accuracies_a)) + 1;



% Coherence
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification_across1.mat';
load([results_directory results_filename]);
coherencies = permute(mean(coherencies(freq_range_a, :, :, :, :), 1), [5 3 4 1 2]); % mean across frequency range; flies x ch-pairs x conditions
scores = zeros(size(coherencies, 2), 1); % ch-pairs
for pair = 1 : size(coherencies, 2)
    % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
    [h, p, ci, stats] = ttest(coherencies(:, pair, 1), coherencies(:, pair, 2));
    [p, h, stats] = signrank(coherencies(:, pair, 1), coherencies(:, pair, 2));
    scores(pair) = stats.signedrank; % stats.signedrank; stats.tstat;
end
% Add to boxplot data structure
measure_accuracies_a = [measure_accuracies_a; scores];
measure_groups_a = [measure_groups_a; zeros(size(scores)) + 2];



% Phi-three
[phis, measure_string] = phi_load('phi_three', 1, bin_location);
group_counter = 3;
for nChannels_counter = 1 : length(phis)
    values = permute(phis{nChannels_counter}.phis(:, :, :, :, tau), [3 1 4 2 5]); % flies x sets x conditions
    phis{nChannels_counter}.tscores = zeros(size(values, 2), 1); % ch-sets
    for set_counter = 1 : size(values, 2)
            % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
            [h, p, ci, stats] = ttest(values(:, set_counter, 1), values(:, set_counter, 2));
            [p, h, stats] = signrank(values(:, set_counter, 1), values(:, set_counter, 2));
            phis{nChannels_counter}.tscores(set_counter) = stats.signedrank; %stats.signedrank; stats.tstat;
    end
    % Add to boxplot data structure
    scores = phis{nChannels_counter}.tscores;
    measure_accuracies_a = [measure_accuracies_a; scores];
    measure_groups_a = [measure_groups_a; zeros(size(scores)) + group_counter];
    group_counter = group_counter + 1;
end
% Store for other plots (store in accuracies - in order to reuse code)
if strcmp(measure, 'phi_three')
    accuracies_a = cell(size(phis));
    for nChannels_counter = 1 : length(phis)
        accuracies_a{nChannels_counter}.accuracies = phis{nChannels_counter}.tscores;
        accuracies_a{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
    end
end



% Phi-star
[phis, measure_string] = phi_load('phi_star', 1, bin_location);
for nChannels_counter = 1 : length(phis)
    values = permute(phis{nChannels_counter}.phis(:, :, :, :, tau), [3 1 4 2 5]); % flies x sets x conditions
    phis{nChannels_counter}.tscores = zeros(size(values, 2), 1); % ch-sets
    for set_counter = 1 : size(values, 2)
            % ttest(a, b) = ttest(a-b), so positive t-value indicates a>b
            [h, p, ci, stats] = ttest(values(:, set_counter, 1), values(:, set_counter, 2));
            [p, h, stats] = signrank(values(:, set_counter, 1), values(:, set_counter, 2));
            phis{nChannels_counter}.tscores(set_counter) = stats.signedrank; %stats.signedrank; stats.tstat;
    end
    % Add to boxplot data structure
    scores = phis{nChannels_counter}.tscores;
    measure_accuracies_a = [measure_accuracies_a; scores];
    measure_groups_a = [measure_groups_a; zeros(size(scores)) + group_counter];
    group_counter = group_counter + 1;
end
% Store for other plots (store in accuracies - in order to reuse code)
if strcmp(measure, 'phi_star')
    accuracies_a = cell(size(phis));
    for nChannels_counter = 1 : length(phis)
        accuracies_a{nChannels_counter}.accuracies = phis{nChannels_counter}.tscores;
        accuracies_a{nChannels_counter}.channel_sets = phis{nChannels_counter}.channel_sets;
    end
end

%% Save

save(out_file,...
    'measure_accuracies_w', 'measure_groups_w',...
    'measure_accuracies_a', 'measure_groups_a',...
    'accuracies_w', 'accuracies_a'...
    );
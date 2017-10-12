%% DESCRIPTION

%{
This script tests the second hypothesis: phi-star is correlated with phi-three

%}

%% SETUP

star_metric = 'phi_stars';

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

filter_percent = 5; % 0 to 100 %

condition_shapes{1} = 'o'; condition_shapes{2} = 'x';
tau_colours = {'r', 'g', 'b'};
tau_alphas = [1 0.8 0.5];

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    %'_shareFiltered'
    ];

share_pairs = 0;

%% LOAD

% disp('loading');
% 
% % Phi-3
% load([data_directory data_filename '_phithree.mat']);
% phi_threes = phis;
% 
% % Phi-star
% load([data_directory data_filename '_phistar.mat']);
% phi_stars = phis;
% if numel(phi_stars) == length(phi_stars) && size(phi_stars, 1) ~= 1
%     phi_stars = phi_stars';
% end
% 
% disp('loaded');

%% Filter for non-sharing/sharing channel sets

% First row contains sets with no pairs of channels sharing an electrode
% Second row contains sets with at least one pair of channels sharing an electrode
% So use option+1: (0+1)=no-share, (1+1)=share

phi_threes = phi_threes(share_pairs+1, :);
phi_stars = phi_stars(share_pairs+1, :);

%% Extract top and bottom 5% of phi values, and top/bottom 5% of delta phi
% % We sort both phi-3 and phi-star using sorted phi-3 indexes
% % For delta phi, we want the phis of the channel sets which have greatest/lowest delta
% 
% nChannels = length(phi_threes);
% 
% for nChannels_counter = 1 : nChannels
%     
%     % Determine range to take based on filter_percent
%     range_end = ceil(filter_percent * (size(phi_threes{nChannels_counter}.phi_threes, 1) / 100));
%     
%     % Sort either phi-three or phi-star, then use sorting index to extract related phi values of the other
%     
%     % Top
%     [phi_threes{nChannels_counter}.top_phis, indices] = sort_linear_index(phi_threes{nChannels_counter}.phi_threes, 1, 'descend');
%     phi_stars{nChannels_counter}.top_phis = phi_stars{nChannels_counter}.phi_stars(indices);
%     phi_threes{nChannels_counter}.indices_descend = indices;
%     phi_stars{nChannels_counter}.indicies_descend = indices;
%     
%     phi_threes{nChannels_counter}.top_phis = phi_threes{nChannels_counter}.top_phis(1:range_end, :, :, :, :);
%     phi_stars{nChannels_counter}.top_phis = phi_stars{nChannels_counter}.top_phis(1:range_end, :, :, :, :);
%     
%     % Bottom
%     [phi_threes{nChannels_counter}.bot_phis, indices] = sort_linear_index(phi_threes{nChannels_counter}.phi_threes, 1, 'ascend');
%     phi_stars{nChannels_counter}.bot_phis = phi_stars{nChannels_counter}.phi_stars(indices);
%     phi_threes{nChannels_counter}.indices_ascend = indices;
%     phi_stars{nChannels_counter}.indicies_ascend = indices;
%     
%     phi_threes{nChannels_counter}.bot_phis = phi_threes{nChannels_counter}.bot_phis(1:range_end, :, :, :, :);
%     phi_stars{nChannels_counter}.bot_phis = phi_stars{nChannels_counter}.bot_phis(1:range_end, :, :, :, :);
%     
%     % Delta
%     phi_threes{nChannels_counter}.delta_phis = squeeze(phi_threes{nChannels_counter}.phi_threes(:, :, :, 1, :) - phi_threes{nChannels_counter}.phi_threes(:, :, :, 2, :));
%     phi_stars{nChannels_counter}.delta_phis = squeeze(phi_stars{nChannels_counter}.phi_stars(:, :, :, 1, :) - phi_stars{nChannels_counter}.phi_stars(:, :, :, 2, :));
%     
%     % Top
%     [phi_threes{nChannels_counter}.top_deltas, indices] = sort_linear_index(phi_threes{nChannels_counter}.delta_phis, 1, 'descend');
%     indices = permute(indices, [1 2 3 5 4]); % Re add condition dimension
%     indices = cat(4, indices, indices);
%     phi_threes{nChannels_counter}.top_deltas = phi_threes{nChannels_counter}.phi_threes();
%     phi_stars{nChannels_counter}.top_deltas = phi_stars{nChannels_counter}.phi_stars(indices);
%     phi_threes{nChannels_counter}.delta_indices_descend = indices;
%     phi_stars{nChannels_counter}.delta_indicies_descend = indices;
%     
%     phi_threes{nChannels_counter}.top_deltas = phi_threes{nChannels_counter}.top_deltas(1:range_end, :, :, :, :);
%     phi_stars{nChannels_counter}.top_deltas = phi_stars{nChannels_counter}.top_deltas(1:range_end, :, :, :, :);
%     
%     % Bottom
%     [phi_threes{nChannels_counter}.bot_deltas, indices] = sort_linear_index(phi_threes{nChannels_counter}.delta_phis, 1, 'ascend');
%     indices = permute(indices, [1 2 3 5 4]); % Re add condition dimension
%     indices = cat(4, indices, indices);
%     phi_threes{nChannels_counter}.bot_deltas = phi_threes{nChannels_counter}.phi_threes(indices);
%     phi_stars{nChannels_counter}.bot_deltas = phi_stars{nChannels_counter}.phi_stars(indices);
%     phi_threes{nChannels_counter}.delta_indices_ascend = indices;
%     phi_stars{nChannels_counter}.delta_indicies_ascend = indices;
%     
%     phi_threes{nChannels_counter}.bot_deltas = phi_threes{nChannels_counter}.bot_deltas(1:range_end, :, :, :, :);
%     phi_stars{nChannels_counter}.bot_deltas = phi_stars{nChannels_counter}.bot_deltas(1:range_end, :, :, :, :);
% end
% 
% %% Average across channel sets
% 
% % for nChannels_counter = 1 : nChannels
% %     phi_threes{nChannels_counter}.phi_threes_setMean = mean(phi_threes{nChannels_counter}.phi_threes, 1);
% %     phi_stars{nChannels_counter}.phi_stars_setMean = mean(phi_stars{nChannels_counter}.phi_stars, 1);
% % end
% 
% %% Correlate for each nChannels and tau lag
% 
% % % Correlation considers every channel combination and trial
% % 
% % nChannels = length(phi_threes); % phi_threes will always have fewer nChannels than phi_stars, due to computability
% % 
% % correlation_results = cell(nChannels, 1);
% % 
% % % One subplot per nChannels
% % % For each subplot: 3 colours (one per tau lag)
% % figure;
% % for nChannels_counter = 1 : nChannels
% %     % Assumes that phi_threes{x}.nChannels is always the same as phi_stars{x}.nChannels
% %     channels_used = phi_threes{nChannels_counter}.nChannels;
% %     taus = phi_threes{nChannels_counter}.taus;
% %     
% %     correlation_results{nChannels_counter}.nChannels = channels_used;
% %     correlation_results{nChannels_counter}.taus = taus;
% %     correlation_results{nChannels_counter}.correlation = zeros(1, length(taus));
% %     correlation_results{nChannels_counter}.correlation_ps = zeros(1, length(taus));
% %     
% %     subplot(1, nChannels, nChannels_counter);
% %     plot_flies = [1 2 3 4 5 6 7 8 9 10 12 13];
% %     
% %     % Correlate for each tau
% %     for tau_counter = 1 : length(taus)
% % %         phi_three_values = log(phi_threes{nChannels_counter}.phi_threes(:, :, :, :, tau_counter));
% % %         phi_star_values = log(phi_stars{nChannels_counter}.phi_stars(:, :, :, :, tau_counter));
% %         phi_three_values = (phi_threes{nChannels_counter}.top_deltas(:, :, 13, 1, tau_counter));
% %         phi_star_values = (phi_stars{nChannels_counter}.top_deltas(:, :, 13, 1, tau_counter));
% %         
% %         %phi_three_values = (phi_threes{nChannels_counter}.top_phis(:, :, 1, 1, tau_counter));
% %         %phi_star_values = (phi_threes{nChannels_counter}.top_phis(:, :, 1, 2, tau_counter));
% %         
% %         % Correlate
% %         [correlation, p] = corrcoef(phi_three_values(:), phi_star_values(:));
% %         % Save (always comparing 2 things, so correlation and p are always in position (1, 2) (or 2, 1)
% %         correlation_results{nChannels_counter}.correlation(tau_counter) = correlation(1, 2);
% %         correlation_results{nChannels_counter}.correlation_ps(tau_counter) = p(1, 2);
% %         
% %         % Plot
% %         scatter(phi_three_values(:), phi_star_values(:), 100, ['.' tau_colours{tau_counter}], 'MarkerEdgeAlpha', tau_alphas(tau_counter));
% %         %axis([-24 0 -24 0]);
% %         %axis([0 0.1 0 0.00002]);
% %         
% %         hold on;
% %     end
% % end
% 
%% Correlation after averaging values across trials

% % Normalise
% for nChannels_counter = 1 : length(phi_threes)
%     phi_threes{nChannels_counter}.phi_threes = zscore(phi_threes{nChannels_counter}.phi_threes, [], 3);
%     phi_stars{nChannels_counter}.phi_stars = zscore(phi_stars{nChannels_counter}.phi_stars, [], 3);
% end

% Average values across trials
phi_threes_avg = cell(length(phi_threes), 1);
phi_stars_avg = cell(length(phi_stars), 1);
for nChannels_counter = 1 : length(phi_threes)
    phi_threes_avg{nChannels_counter} = (squeeze(mean(phi_threes{nChannels_counter}.phi_threes, 2)));
    phi_stars_avg{nChannels_counter} = (squeeze(mean(phi_stars{nChannels_counter}.(star_metric), 2)));
    %phi_stars_avg{nChannels_counter} = (squeeze(mean(phi_stars{nChannels_counter}.mis - phi_stars{nChannels_counter}.mi_stars, 2)));
end

% Correlate for each fly, condition, tau level
correlations = cell(length(phi_threes), 1);
for nChannels_counter = 1 : length(phi_threes_avg)
    correlations{nChannels_counter} = zeros(size(phi_threes_avg{nChannels_counter}, 2), size(phi_threes_avg{nChannels_counter}, 3), size(phi_threes_avg{nChannels_counter}, 4));
end
figure;
subplot_counter = 1;
for tau_counter = 1 : size(phi_threes_avg{nChannels_counter}, 4)
    for nChannels_counter = 1 : length(phi_threes_avg)
        subplot(length(phi_threes{1}.taus), length(phi_threes_avg), subplot_counter);
        for fly_counter = 1 : size(phi_threes_avg{nChannels_counter}, 2)
            for condition = 1 : size(phi_threes_avg{nChannels_counter}, 3)
                three = phi_threes_avg{nChannels_counter}(:, fly_counter, condition, tau_counter);
                star = phi_stars_avg{nChannels_counter}(:, fly_counter, condition, tau_counter);
                
                %scatter(three(:), star(:), 100, [condition_shapes{condition} tau_colours{tau_counter}], 'MarkerEdgeAlpha', tau_alphas(tau_counter));
                %axis([0 0.25 0 0.018]);
                %hold on;
                
                scatter(three(:), star(:), condition_shapes{condition}); hold on;
                set(gca, 'xscale', 'log');
                set(gca, 'yscale', 'log'); % correlation test will be conducted after log transform
                axis tight
                
                if subplot_counter <= length(phi_threes_avg)
                    title([num2str(phi_threes{nChannels_counter}.nChannels) ' channels']);
                end
                if subplot_counter == (length(phi_threes{1}.taus) * length(phi_threes_avg)) - (length(phi_threes_avg)-1)
                    xlabel('phi-3');
                    ylabel(['phi-* tau=' num2str(phi_threes{1}.taus(tau_counter))]);
                elseif mod(subplot_counter, length(phi_threes_avg)) == 1
                    ylabel(['tau=' num2str(phi_threes{1}.taus(tau_counter))]);
                end
                
                % Compute and store correlation
                [correlation, p] = corrcoef((three), (star));
                correlations{nChannels_counter}(fly_counter, condition, tau_counter) = correlation(1, 2);
            end
        end
        subplot_counter = subplot_counter + 1;
    end
    
end

% Plot correlations per fly
figure;
for nChannels_counter = 1 : length(correlations)
    subplot(1, 3, nChannels_counter);
    
    for fly_counter = 1 : size(correlations{nChannels_counter}, 1)
        for condition_counter = 1 : size(correlations{nChannels_counter}, 2)
            for tau_counter = 1 : size(correlations{nChannels_counter}, 3)
                scatter(fly_counter, correlations{nChannels_counter}(fly_counter, condition_counter, tau_counter), [condition_shapes{condition_counter} tau_colours{tau_counter}]);
                axis([0 size(correlations{nChannels_counter}, 1)+1 0 1]);
                hold on;
                
                title([num2str(phi_threes{nChannels_counter}.nChannels) ' channels']);
                if nChannels_counter == 1
                    xlabel('fly');
                    ylabel('r');
                end
            end
        end
    end

end

%% Plot figure
% Top row: relationship between phis for one fly
% Bottom row: coeffs for all flies
plot_fly = 1;
condition_colours = 'rb';
plot_tau = 1;
tau_shapes = '*o^';

figure;
for nChannels_counter = 1 : length(correlations)
    subplot(2, length(correlations), nChannels_counter);
    for tau_counter = 1 : plot_tau %size(phi_threes_avg{nChannels_counter}, 4)
        for condition_counter = 1 : size(phi_threes_avg{nChannels_counter}, 3)
            three = phi_threes_avg{nChannels_counter}(:, plot_fly, condition_counter, tau_counter);
            star = phi_stars_avg{nChannels_counter}(:, plot_fly, condition_counter, tau_counter);
            
            scatter(three(:), star(:), [], [condition_shapes{condition_counter}]); hold on;
        end
        set(gca, 'xscale', 'log');
        set(gca, 'yscale', 'log'); % correlation test will be conducted after log transform
        axis tight
        title([num2str(phi_threes{nChannels_counter}.nChannels) ' channels']);
        ylabel([char(981) '*'], 'FontSize', 15, 'rotation', 0);
        xlabel(char(981), 'FontSize', 15, 'rotation', 0);
    end
end

% Correlations for each fly
correlations = cell(length(phi_threes), 1);
correlations_ps = cell(length(phi_threes), 1);
for nChannels_counter = 1 : length(phi_threes_avg)
    correlations{nChannels_counter} = zeros(size(phi_threes_avg{nChannels_counter}, 2), size(phi_threes_avg{nChannels_counter}, 4));
    correlations_ps{nChannels_counter} = zeros(size(phi_threes_avg{nChannels_counter}, 2), size(phi_threes_avg{nChannels_counter}, 4));
end
for nChannels_counter = 1 : length(correlations)
    for fly_counter = 1 : size(phi_threes_avg{nChannels_counter}, 2)
        for tau_counter = 1 : size(phi_threes_avg{nChannels_counter}, 4)
            three = phi_threes_avg{nChannels_counter}(:, fly_counter, :, tau_counter);
            star = phi_stars_avg{nChannels_counter}(:, fly_counter, :, tau_counter);
            [correlation, p] = corrcoef((three(:)), (star(:)));
            correlations{nChannels_counter}(fly_counter, tau_counter) = correlation(1, 2);
            correlations_ps{nChannels_counter}(fly_counter, tau_counter) = p(1, 2);
        end
    end
end

for nChannels_counter = 1 : length(correlations)
    subplot(2, length(correlations), nChannels_counter + length(correlations));
    
    for fly_counter = 1 : size(correlations{nChannels_counter}, 1)
        for tau_counter = 1 : size(correlations{nChannels_counter}, 2)
            scatter(fly_counter, correlations{nChannels_counter}(fly_counter, tau_counter), [tau_shapes(tau_counter) tau_colours{tau_counter}]);
            axis([0 size(correlations{nChannels_counter}, 1)+1 0 1]);
            hold on;
            
            if nChannels_counter == 1
                xlabel('fly', 'FontSize', 15, 'rotation', 0);
                ylabel('r', 'FontSize', 15, 'rotation', 0);
            end
        end
    end
    
end

%% Correlations after matching MIPs

%mipmatch_correlation(phi_threes, phi_threes_avg, phi_stars_avg, condition_shapes, 1, tau_colours);

%% Correlations per set

[correlation_rs, correlations_ps, correlation_sigs] = correlate_per_set(phi_threes_avg, phi_stars_avg);

%% Function: correlation after filtering for 'matching' MIPs

function [] = mipmatch_correlation(phi_threes, phi_threes_avg, phi_stars_avg, condition_shapes, plot_fly, tau_colours)
% Plots correlations after filtering for matching MIPs
% Mips 'match' if the portion of matching is at least 2*chance
%
% Inputs:
%   phi_threes: phi_threes, averaged across trials
%   phi_stars: phi_stars, averaged across trials

chance_levels = [1 1/3 1/7];
chance_multiplier_threshold = 2;
samples_per_trial = 2250;

% Load match results
data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;
results_directory = 'analysis_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_mipComparison.mat'
    ];
disp('loading matches');
load([results_directory results_filename]);
disp('loaded');

% Find sets which meet the above chance threshold
match_sets = cell(length(matches_per_trial), 1);
for nChannels_counter = 1 : length(matches_per_trial)
    dims = size(matches_per_trial{nChannels_counter}.matches);
    match_sets{nChannels_counter} = zeros([dims(1) dims(3:end)]);
    for tau = 1 : size(matches_per_trial{nChannels_counter}.matches, 5)
        for condition = 1 : size(matches_per_trial{nChannels_counter}.matches, 4)
            for fly = 1 : size(matches_per_trial{nChannels_counter}.matches, 3)
                for channel_set = 1 : size(matches_per_trial{nChannels_counter}.matches, 1)
                    trial_matches = sum(matches_per_trial{nChannels_counter}.matches(channel_set, :, fly, condition, tau));
                    if trial_matches >= samples_per_trial*size(matches_per_trial{nChannels_counter}.matches, 2)*(chance_levels(nChannels_counter) * chance_multiplier_threshold)
                        match_sets{nChannels_counter}(channel_set, fly, condition, tau) = 1;
                    end
                end
            end
        end
    end
end

% Find correlations per fly
% We will make the plots at all parameters for 1 fly
% We will calculate the correlation for all flies
figure;
subplot_counter = 1;
for tau_counter = 1 : size(phi_threes_avg{nChannels_counter}, 4)
    for nChannels_counter = 2 : length(phi_threes_avg)
        subplot(length(phi_threes{1}.taus), length(phi_threes_avg)-1, subplot_counter);
        for fly_counter = plot_fly : plot_fly
            for condition = 1 : size(phi_threes_avg{nChannels_counter}, 3)
                three = phi_threes_avg{nChannels_counter}(not(logical(match_sets{nChannels_counter}(:, fly_counter, condition, tau_counter))), fly_counter, condition, tau_counter);
                star = phi_stars_avg{nChannels_counter}(not(logical(match_sets{nChannels_counter}(:, fly_counter, condition, tau_counter))), fly_counter, condition, tau_counter);
                
                scatter(three(:), star(:), condition_shapes{condition}); hold on;
                set(gca, 'xscale', 'log');
                set(gca, 'yscale', 'log');
                
                if subplot_counter <= length(phi_threes_avg)
                    title([num2str(phi_threes{nChannels_counter}.nChannels) ' channels']);
                end
                if subplot_counter == (length(phi_threes{1}.taus) * (length(phi_threes_avg)-1)) - 1
                    xlabel('phi-3');
                    ylabel(['phi-* tau=' num2str(phi_threes{1}.taus(tau_counter))]);
                elseif mod(subplot_counter, length(phi_threes_avg)-1) == 1
                    ylabel(['tau=' num2str(phi_threes{1}.taus(tau_counter))]);
                end
            end
        end
        subplot_counter = subplot_counter + 1;
    end
end


% Calculate correlation for each fly
correlations = cell(length(phi_threes), 1);
for nChannels_counter = 1 : length(phi_threes_avg)
    correlations{nChannels_counter} = zeros(size(phi_threes_avg{nChannels_counter}, 2), size(phi_threes_avg{nChannels_counter}, 3), size(phi_threes_avg{nChannels_counter}, 4));
end
for tau_counter = 1 : size(phi_threes_avg{nChannels_counter}, 4)
    for nChannels_counter = 2 : length(phi_threes_avg)
        for fly_counter = 1:size(phi_threes_avg{nChannels_counter}, 2)
            for condition = 1 : size(phi_threes_avg{nChannels_counter}, 3)
                three = phi_threes_avg{nChannels_counter}(not(logical(match_sets{nChannels_counter}(:, fly_counter, condition, tau_counter))), fly_counter, condition, tau_counter);
                star = phi_stars_avg{nChannels_counter}(not(logical(match_sets{nChannels_counter}(:, fly_counter, condition, tau_counter))), fly_counter, condition, tau_counter);
                % Compute and store correlation
                [correlation, p] = corrcoef(three, star);
                correlations{nChannels_counter}(fly_counter, condition, tau_counter) = correlation(1, 2);
            end
        end
    end
end

% Plot correlations per fly
figure;
for nChannels_counter = 2 : length(correlations)
    subplot(1, 2, nChannels_counter-1);
    
    for fly_counter = 1 : size(correlations{nChannels_counter}, 1)
        for condition_counter = 1 : size(correlations{nChannels_counter}, 2)
            for tau_counter = 1 : size(correlations{nChannels_counter}, 3)
                scatter(fly_counter, correlations{nChannels_counter}(fly_counter, condition_counter, tau_counter), [condition_shapes{condition_counter} tau_colours{tau_counter}]);
                axis([0 size(correlations{nChannels_counter}, 1)+1 0 1]);
                hold on;
                
                title([num2str(phi_threes{nChannels_counter}.nChannels) ' channels']);
                if nChannels_counter == 2
                    xlabel('fly');
                    ylabel('r');
                end
            end
        end
    end

end

end

%% Function: correlate per channel set

function [correlations, correlation_ps, correlation_sigs] = correlate_per_set(phi_threes_avg, phi_stars_avg)
% Takes trial averaged values (i.e. one value per set, fly, condition, tau)

q = 0.05;
height_sig = 1;
height_nonsig = 0.1;
nChannels_list = [2 3 4];
tau_list = [4 8 16];


correlations = cell(length(phi_threes_avg), 1);
correlation_ps = cell(length(phi_threes_avg), 1);
for nChannels_counter = 1 : length(phi_threes_avg)
    correlations{nChannels_counter} = zeros(size(phi_threes_avg{nChannels_counter}, 1), size(phi_threes_avg{nChannels_counter}, 3), size(phi_threes_avg{nChannels_counter}, 4));
    correlation_ps{nChannels_counter} = zeros(size(phi_threes_avg{nChannels_counter}, 1), size(phi_threes_avg{nChannels_counter}, 3), size(phi_threes_avg{nChannels_counter}, 4));
end

% Compute correlations, combining across nChannels
for nChannels_counter = 1 : length(phi_threes_avg)
    for tau = 1 : size(phi_threes_avg{nChannels_counter}, 4)
        for condition = 1 : size(phi_threes_avg{nChannels_counter}, 3)
            for set = 1 : size(phi_threes_avg{nChannels_counter}, 1)
                three = phi_threes_avg{nChannels_counter}(set, :, condition, tau);
                star = phi_stars_avg{nChannels_counter}(set, :, condition, tau);
                [correlation, p] = corrcoef((three), (star));
                correlations{nChannels_counter}(set, condition, tau) = correlation(1, 2);
                correlation_ps{nChannels_counter}(set, condition, tau) = p(1, 2);
            end
        end
    end
end

% Concatenate nChannels (to get a single axis)
correlations_all = [];
correlation_ps_all = [];
channel_ticks = []; xtick_counter = 1;
channel_labels = [];
for nChannels_counter = 1 : length(correlations)
    correlations_all = cat(1, correlations_all, correlations{nChannels_counter});
    correlation_ps_all = cat(1, correlation_ps_all, correlation_ps{nChannels_counter});
    channel_ticks = cat(1, channel_ticks, xtick_counter);
    channel_labels = cat(1, channel_labels, nChannels_list(nChannels_counter));
    xtick_counter = xtick_counter + size(correlations{nChannels_counter}, 1);
end

% Plot
figure;
subplot_counter = 1;
correlation_sigs = zeros(size(correlation_ps_all, 1), size(correlations_all, 2), size(correlations_all, 3));
for tau = 1 : size(correlations_all, 3)
    for condition = 1 : size(correlations_all, 2)
        subplot(size(correlations_all, 3), size(correlations_all, 2), subplot_counter);
        
        bar(correlations_all(:, condition, tau));
        axis([-50 size(correlations_all, 1)+50 -0.3 1.1]);
        hold on;
        
        % Plot significance after FDR correction
        corrected = fdr_correct(correlation_ps_all(:, condition, tau), q);
        xs = (1:length(corrected));
        scatter(xs(logical(corrected)), corrected(logical(corrected))*height_sig, 10, 'k*', 'MarkerEdgeAlpha', 0.1); % Sig
        scatter(xs(not(logical(corrected))), (corrected(not(logical(corrected))))-height_nonsig, 10, 'k*', 'MarkerEdgeAlpha', 0.1); % Not sig
        correlation_sigs(:, condition, tau) = corrected; % Store for returning
        
        xticks(channel_ticks); xticklabels(channel_labels);
        
        if subplot_counter == 1
            title('air');
        elseif subplot_counter == 2
            title('iso');
        end
        
        if subplot_counter == (size(correlations_all, 3) * size(correlations_all, 2)) - 1
            ylabel(['r tau=' num2str(tau_list(tau))]);
        elseif mod(subplot_counter, size(correlations_all, 2)) == 1
            ylabel(['tau=' num2str(tau_list(tau))]);
        end
        
        subplot_counter = subplot_counter + 1;
    end
end

end

%% DESCRIPTION

%{
This script takes the relevant extracts data corresponding to rest periods
(i.e. periods which have an air-puff)
%}

%% SETUP

prep_detrend = 0;
prep_zscore = 0;

fly = 1;
channels = [2 4];

taus = [1 2 4 8 16];
epoch_splits= [1 2 4 8 16];
% 2.25s epoch splits (truncate where needed to accommodate tau):
% 1 = 2250ms
% 2 = 1125
% 4 = 562.5
% 8 = 281.25
% 16 = 140.625

data_directory = 'workspace_results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim.mat';

results_directory = 'workspace_results/';
results_filename = [data_file '_detrend' num2str(prep_detrend) '_zscore' num2str(prep_zscore) '_fly' num2str(fly) 'chs' num2str(channels(1)) 'and' num2str(channels(2))];

condition_colours = {'r', 'b'};
condition_markers = {'o', 'x'};

%% LOAD

load([data_directory data_file]);

%% START

%data = permute(fly_data(:, channels, :, fly, :), [1 2 3 5 4]); % samples x channels x trials x conditions

if prep_detrend == 1
    for condition = 1 : 2
        for trial = 1 : size(fly_data, 3)
            fly_data(:, :, trial, fly, condition) = detrend(fly_data(:, :, trial, fly, condition));
        end
    end
end
if prep_zscore == 1
    fly_data = zscore(fly_data, [], 1);
end

results = cell(1, 2); % 1=air, 2=anest
results{1} = cell(length(epoch_splits), length(taus));
results{2} = cell(length(epoch_splits), length(taus));

for split_counter = 1 : length(epoch_splits)
    nBins = epoch_splits(split_counter);
    for tau_counter = 1 : length(taus)
        tau = taus(tau_counter);
        
        for condition = 1 : 2
            results{condition}{split_counter, tau_counter}.nBins = nBins;
            results{condition}{split_counter, tau_counter}.tau = tau;
            results{condition}{split_counter, tau_counter}.cov_present_present_det = zeros(nBins*size(fly_data, 3), 1);
            results{condition}{split_counter, tau_counter}.cov_present_past_det = zeros(nBins*size(fly_data, 3), 1);
            results{condition}{split_counter, tau_counter}.cov_past_past_det = zeros(nBins*size(fly_data, 3), 1);
            results{condition}{split_counter, tau_counter}.phistar = zeros(nBins*size(fly_data, 3), 2);
%             results{condition}{split_counter, tau_counter}.condh = zeros(nBins*size(fly_data, 3), 2); % If using Masafumi's phi_gauss.m
%             results{condition}{split_counter, tau_counter}.mi = zeros(nBins*size(fly_data, 3), 2); % If using Masafumi's phi_gauss.m
            results{condition}{split_counter, tau_counter}.condh = zeros(nBins*size(fly_data, 3), 1);
            results{condition}{split_counter, tau_counter}.mi = zeros(nBins*size(fly_data, 3), 1);
            results{condition}{split_counter, tau_counter}.h = zeros(nBins*size(fly_data, 3), 1);
            
            trialBin_counter = 1;
            for trial = 1 : size(fly_data, 3)
                
                % Split epoch into nBins
                data_bins = split_data(fly_data(:, channels, trial, fly, condition), nBins, tau);
                
                % Calculate phi per bin
                for bin_counter = 1 : nBins
                    bin = [data_bins(1:tau, :, 1, bin_counter)' data_bins(:, :, 2, bin_counter)'];
                    [cov_present_present, cov_present_past, cov_past_past] = Cov_comp_sample(bin, tau);
                    %[cov_present_present, cov_present_past, ~, cov_past_past] = Cov_comp_shrink(bin, tau);
                    
                    [...
                        results{condition}{split_counter, tau_counter}.phistar(trialBin_counter, :),...
                        ~, ~, ~,...
                        results{condition}{split_counter, tau_counter}.mi(trialBin_counter),...
                        ~, ~,...
                        results{condition}{split_counter, tau_counter}.condh(trialBin_counter),...
                        results{condition}{split_counter, tau_counter}.h(trialBin_counter)] = phi_comp(cov_present_present, cov_present_past, cov_past_past);
                    
%                     [...
%                         results{condition}{split_counter, tau_counter}.phistar(trialBin_counter, :),...
%                         results{condition}{split_counter, tau_counter}.mi(trialBin_counter, :),...
%                         results{condition}{split_counter, tau_counter}.condh(trialBin_counter, :)...
%                         ] = phi_gauss(cov_present_present, cov_present_past, [], cov_past_past, []);
                    
                    results{condition}{split_counter, tau_counter}.cov_present_present_det(trialBin_counter) = det(cov_present_present);
                    results{condition}{split_counter, tau_counter}.cov_present_past_det(trialBin_counter) = det(cov_present_past);
                    results{condition}{split_counter, tau_counter}.cov_past_past_det(trialBin_counter) = det(cov_past_past);
                    
                    trialBin_counter = trialBin_counter + 1;
                    
                end
                
            end
        end
    end
end

%% Plot (mean + std)

% cov_present_present_det
% cov_present_past_det
% cov_past_past_det
% phistar (note you need to plot (:, x); x=1 is non-normalised, x=2 is)
% condh (conditional entropy)
% h (entropy)
% mi (mutual information)
plot_metric = 'phistar';

result_stats = struct();
result_stats.mean = zeros(size(results{1}, 1), 2, length(taus));
result_stats.std = zeros(size(results{1}, 1), 2, length(taus));
result_stats.var = zeros(size(results{1}, 1), 2, length(taus));

% Summary stats
for tau_counter = 1 : length(taus)
    
    for split_counter = 1 : length(epoch_splits)
        
        for condition = 1 : 2
            result_stats.mean(split_counter, condition, tau_counter) = mean(results{condition}{split_counter, tau_counter}.(plot_metric)(:, 1));
            result_stats.std(split_counter, condition, tau_counter) = std(results{condition}{split_counter, tau_counter}.(plot_metric)(:, 1));
            result_stats.var(split_counter, condition, tau_counter) = var(results{condition}{split_counter, tau_counter}.(plot_metric)(:, 1));
        end
    end
    
end

% Plot
figure;
condition_xs = zeros(length(epoch_splits), 2);
condition_xs(:, 1) = (1:length(epoch_splits)) - 0.1;
condition_xs(:, 2) = (1:length(epoch_splits)) + 0.1;
for tau_counter = 1 : length(taus)
    subplot(2, 3, tau_counter);
    
    for condition = 1 : 2
        errorbar(condition_xs(:, condition), result_stats.mean(:, condition), result_stats.std(:, condition, tau_counter), condition_colours{condition}); hold on;
    end
    
    title(['tau = ' num2str(taus(tau_counter))]);
    set(gca, 'XTick', [1, 2, 3, 4, 5], 'XTickLabel', [1, 2, 4, 8, 16]);
    
    if tau_counter == 1
        ylabel([plot_metric]);
    elseif tau_counter == 5
        xlabel('bins');
    end
    
    axis tight;
    
    % No shrink axis
    %axis([0 6 -0.01 0.09]); % covXY
    %axis([0 6 0.85 1.05]); % covX
    %axis([0 6 -0.01 0.12]); % mi
    %axis([0 6 2.75 2.86]); % h
    %axis([0 6 -5*10^-3 18*10^-3]); % phistar
    
    % Shrink axis
    %axis([0 6 -1*10^4 3*10^4]); % covXY
    %axis([0 6 1.4*10^5 3.1*10^5]); % covX
    %axis([0 6 -0.01 0.12]); % mi
    %axis([0 6 8.8 9.2]); % h
    axis([0 6 -5*10^-3 18*10^-3]); % phistar
end

%%
%print([plot_metric '_shrink'], '-dpng');

%% Plot (raw)

% cov_present_present_det
% cov_present_past_det
% cov_past_past_det
% phistar (note you need to plot (:, x); x=1 is non-normalised, x=2 is)
% condh (conditional entropy)
% mi (mutual information)
% plot_metric = 'cov_present_past_det';
% marker_size = 20;
% 
% figure;
% 
% for tau_counter = 1 : length(taus)
%     subplot(2, 3, tau_counter);
%     
%     for split_counter = 1 : length(epoch_splits)
%         nBins = epoch_splits(split_counter);
%         
%         for condition = 1 : 2
%             scatter(zeros(nBins*size(fly_data, 3), 1)+split_counter, (results{condition}{split_counter, tau_counter}.(plot_metric)(:, 1)), marker_size, condition_colours{condition}, condition_markers{condition}); hold on;
%         end
%     end
%     
%     axis([0 5 -2 0]);
%     title(['tau = ' num2str(taus(tau_counter))]);
%     set(gca, 'XTick', [1, 2, 3, 4, 5], 'XTickLabel', [1, 2, 4, 8, 16]);
%     
%     if tau_counter == 1
%         ylabel([plot_metric]);
%     elseif tau_counter == 5
%         xlabel('bins');
%     end
%     
% end

%%
%print(plot_metric, '-dpng');

%% FUNCTIONS

function [split_bins] = split_data(trial, nBins, tau)
% Splits a trial (samples by trials) into nBin bin pairs (present-past)
% with a sample lag of tau
%
% Inputs:
%   trial = matrix: samples x channels
%   nBins = number of bins to split the trial into
%   tau = sample lag between each present-past bin pair
%
% Outputs:
%   split_data = matrix: samples x channels x 2 (past=1 or present=2) x number
%   of bin pairs;

bin_length = floor(size(trial, 1) / nBins);
split_bins = zeros(bin_length-tau, size(trial, 2), 2, nBins);

bin_position = 1;
for bin_pair = 1 : nBins
    bin = trial(bin_position : bin_position + bin_length - 1, :);
    split_bins(:, :, 1, bin_pair) = bin(1 : end-tau, :);
    split_bins(:, :, 2, bin_pair) = bin(tau+1 : end, :);
    bin_position = bin_position + bin_length;
end

end
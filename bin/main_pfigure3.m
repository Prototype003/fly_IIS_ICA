%% DESCRIPTION

%{

Paper figure 3

Based off of main_phithreevstar_correlations.m and
main_phithreevstar_mipCorrelation_stats.m

%}

%% SETUP

star_metric = 'phi_stars';

taus = [4 8 16];

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

disp('loading');

% Phi-3
load([data_directory data_filename '_phithree.mat']);
phi_threes = phis;

% Phi-star
load([data_directory data_filename '_phistar.mat']);
phi_stars = phis;
if numel(phi_stars) == length(phi_stars) && size(phi_stars, 1) ~= 1
    phi_stars = phi_stars';
end

disp('loaded');

%% Average values across trials

% Average values across trials
phi_threes_avg = cell(length(phi_threes), 1);
phi_stars_avg = cell(length(phi_stars), 1);
for nChannels_counter = 1 : length(phi_threes)
    phi_threes_avg{nChannels_counter} = (squeeze(mean(phi_threes{nChannels_counter}.phi_threes, 2)));
    phi_stars_avg{nChannels_counter} = (squeeze(mean(phi_stars{nChannels_counter}.(star_metric), 2)));
end

%% Calculate correlations per fly

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
            three = phi_threes_avg{nChannels_counter}(:, fly_counter, :, tau_counter); % Correlate across both condition
            star = phi_stars_avg{nChannels_counter}(:, fly_counter, :, tau_counter);
            [correlation, p] = corrcoef((three(:)), (star(:)));
            correlations{nChannels_counter}(fly_counter, tau_counter) = correlation(1, 2);
            correlations_ps{nChannels_counter}(fly_counter, tau_counter) = p(1, 2);
        end
    end
end

%% Average correlations across flies

% Fisher r to z transform
correlations_z = zeros(length(correlations), size(correlations{1}, 2), size(correlations{1}, 1)); % nChannels x taus x flies
correlations_r = zeros(length(correlations), size(correlations{1}, 2), size(correlations{1}, 1)); % original coefficients in new format
for nChannels_counter = 1 : length(correlations)
    for fly_counter = 1 : size(correlations{nChannels_counter}, 1)
        for tau_counter = 1 : size(correlations{nChannels_counter}, 2)
            correlations_z(nChannels_counter, tau_counter, fly_counter) = fisher_rz(correlations{nChannels_counter}(fly_counter, tau_counter));
            correlations_r(nChannels_counter, tau_counter, fly_counter) = correlations{nChannels_counter}(fly_counter, tau_counter);
        end
    end
end

% Mean and standard error
correlations_z_mean = mean(correlations_z, 3);
correlations_z_stderr = std(correlations_z, [], 3) / sqrt(size(correlations_z, 3));

% Fisher z to r transform
correlations_mean = zeros(size(correlations_z_mean));
correlations_stderr = zeros(size(correlations_z_stderr));
for nChannels_counter = 1 : size(correlations_z_mean, 1)
    for tau_counter = 1 : size(correlations_z_mean, 2)
        correlations_mean(nChannels_counter, tau_counter) = fisher_zr(correlations_z_mean(nChannels_counter, tau_counter));
        correlations_stderr(nChannels_counter, tau_counter) = fisher_zr(correlations_z_stderr(nChannels_counter, tau_counter));
    end
end

%% Plot

nChannels_widths = [1 1.5 2];
nChannels_markers = 'x^o';
nChannels_offsets = [-0.5 0 0.5];

figure;
subplot(1, 2, 1);
for nChannels_counter = 1 : size(correlations_mean, 1)
    errorbar(taus+nChannels_offsets(nChannels_counter),...
        correlations_mean(nChannels_counter, :),...
        correlations_stderr(nChannels_counter, :),...
        ['k' nChannels_markers(nChannels_counter) ':'], 'CapSize', 10, 'LineWidth', nChannels_widths(nChannels_counter)); hold on;
end
legend('    2ch', '    3ch', '    4ch', 'Location', 'southeast');
set(gca, 'YTick', [0.4 0.5 0.6], 'FontSize', 12);
set(gca, 'XTick', taus, 'FontSize', 12);
y = ylabel('r', 'FontSize', 15, 'rotation', 90);
xlabel('\tau', 'FontSize', 15);
axis([2 18 0.31 0.64]);

%% Stats

% % ANOVA
% anova_data = []; % flies*nChannels x taus
% for nChannels_counter = 1 : length(correlations)
%     anova_data = [anova_data; permute(correlations_z(nChannels_counter, :, :), [3 2 1])];
% end
% anova_results = struct();
% [anova_results.p, anova_results.table, anova_results.stats] = anova2(anova_data, size(correlations_z, 3));

%% Setup

set_sizes = (3:4);
flies = (1:13);
conditions = (1);
taus = [4 8 16];

conditionSplit = 0;


data_directory = 'analysis_results/';
data_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t4_phithreevstar_correlations_conditionSplit' num2str(conditionSplit) '.mat'];

correlation_matrices = {...
    'phi_v_phistar',...
    'phi_v_mi',...
    'phi_v_mistar',...
    };
p_suffix = '_ps';

%% Load

load([data_directory data_file]);

% Add condition dimension if it's not there already (i.e. conditionSplit is
% 0)
if conditionSplit == 0
    for nChannels = 1 : length(set_sizes)
        for metric_counter = 1 : length(correlation_matrices)
            metric = correlation_matrices{metric_counter};
            correlations{nChannels}.(metric) = permute(correlations{nChannels}.(metric), [1 2 4 3]);
            correlations{nChannels}.([metric p_suffix]) = permute(correlations{nChannels}.([metric p_suffix]), [1 2 4 3]);
            
            % Repeat for correlations on trial averages
            correlations{nChannels}.trial_averaged.(metric) = permute(correlations{nChannels}.trial_averaged.(metric), [1 2 4 3]);
            correlations{nChannels}.trial_averaged.([metric p_suffix]) = permute(correlations{nChannels}.trial_averaged.([metric p_suffix]), [1 2 4 3]);
        end
    end
end

%% Replace non-sig correlations with r = 0

sig_suffix = '_sigs';
q = 0.05;

% Get logical indexes of significant correlations
for nChannels = 1 : length(correlations)
    for metric_counter = 1 : length(correlation_matrices)
        metric = correlation_matrices{metric_counter};
        
        correlations{nChannels}.([metric sig_suffix]) = zeros(size(correlations{nChannels}.(metric)));
        correlations{nChannels}.trial_averaged.([metric sig_suffix]) = zeros(size(correlations{nChannels}.(metric)));
        
        % FDR within each fly and tau
        for fly = 1 : length(flies)
            for condition = 1 : size(correlations{nChannels}.(metric), 3)
                for tau = 1 : length(taus)
                    p_vector = correlations{nChannels}.([metric p_suffix])(:, fly, condition, tau);
                    [threshold_para, threshold_nonpara] = FDR(p_vector, q);
                    if numel(threshold_para) == 0
                        threshold_para = 0;
                    end
                    correlations{nChannels}.([metric sig_suffix])(:, fly, condition, tau) = correlations{nChannels}.([metric p_suffix])(:, fly, condition, tau) < threshold_para;
                    
                    % repeat for trial averaged
                    p_vector = correlations{nChannels}.trial_averaged.([metric p_suffix])(:, fly, condition, tau);
                    [threshold_para, threshold_nonpara] = FDR(p_vector, q);
                    if numel(threshold_para) == 0
                        threshold_para = 0;
                    end
                    correlations{nChannels}.trial_averaged.([metric sig_suffix])(:, fly, condition, tau) = correlations{nChannels}.trial_averaged.([metric p_suffix])(:, fly, condition, tau) < threshold_para;
                end
            end
        end
        
    end
end

% Replace non-significant correlations with 0
for nChannels = 1 : length(correlations)
    for metric_counter = 1 : length(correlation_matrices)
        metric = correlation_matrices{metric_counter};
        
        for element = 1 : numel(correlations{nChannels}.(metric))
            if correlations{nChannels}.([metric sig_suffix])(element) == 0
                correlations{nChannels}.(metric)(element) = 0;
            end
        end
        for element = 1 : numel(correlations{nChannels}.trial_averaged.(metric))
            if correlations{nChannels}.trial_averaged.([metric sig_suffix])(element) == 0
                correlations{nChannels}.trial_averaged.(metric)(element) = 0;
            end
        end
        
    end
end

%% Average correlations across channel sets and flies

% Average ALL (doesn't exclude non-sig) correlations
% Mean all matrices
% Fisher z transform
% Average across channel sets
% Average across flies
% Inverse Fisher z transform

mean_suffix = '_mean';
std_suffix = '_std';

for nChannels = 1 : length(correlations)
    for metric_counter = 1 : length(correlation_matrices)
        metric = correlation_matrices{metric_counter};
        
        % Get correlations
        matrix = correlations{nChannels}.(metric);
        matrix_trial_averaged = correlations{nChannels}.trial_averaged.(metric);
        
        % r to z transform
        for element = 1 : numel(matrix)
            matrix(element) = fisher_rz(matrix(element));
        end
        for element = 1 : numel(matrix_trial_averaged)
            matrix_trial_averaged(element) = fisher_rz(matrix_trial_averaged(element));
        end
        
        % Mean (across sets, then flies)
        correlations{nChannels}.([metric mean_suffix]) = permute(mean(mean(matrix, 1), 2), [4 3 2 1]); % Put the tau dimension first so singleton dimensions are squeezed out
        correlations{nChannels}.trial_averaged.([metric mean_suffix]) = permute(mean(mean(matrix_trial_averaged, 1), 2), [4 3 2 1]);
        
        % Standard error (across flies)
        correlations{nChannels}.([metric std_suffix]) = permute(std(mean(matrix, 1), [], 2), [4 3 2 1]) / sqrt(size(matrix, 2));
        correlations{nChannels}.trial_averaged.([metric std_suffix]) = permute(std(mean(matrix_trial_averaged, 1), [], 2), [4 3 2 1]) / sqrt(size(matrix_trial_averaged, 2));
        
        % z to r transform
        for element = 1 : numel(correlations{nChannels}.([metric mean_suffix]))
            correlations{nChannels}.([metric mean_suffix])(element) = fisher_zr(correlations{nChannels}.([metric mean_suffix])(element));
            correlations{nChannels}.([metric std_suffix])(element) = fisher_zr(correlations{nChannels}.([metric std_suffix])(element));
        end
        for element = 1 : numel(correlations{nChannels}.trial_averaged.([metric std_suffix]))
            correlations{nChannels}.trial_averaged.([metric mean_suffix])(element) = fisher_zr(correlations{nChannels}.trial_averaged.([metric mean_suffix])(element));
            correlations{nChannels}.trial_averaged.([metric std_suffix])(element) = fisher_zr(correlations{nChannels}.trial_averaged.([metric std_suffix])(element));
        end
        
    end
end

%% Plot average correlations

metric = 'phi_v_mistar';
nChannels_colours = {};
nChannels_widths = [1.5 2];
nChannels_markers = '^o';
nChannels_offsets = [-0.4 0.4];

subplot(1, 2, 2);
for nChannels = 1 : length(set_sizes)
    %subplot(1, length(nChannels), nChannels);
    %bar(correlations{nChannels}.trial_averaged.([metric mean_suffix])); hold on;
    
    errorbar(taus+nChannels_offsets(nChannels),...
        correlations{nChannels}.trial_averaged.([metric mean_suffix]),...
        correlations{nChannels}.trial_averaged.([metric std_suffix]),...
        ['k' nChannels_markers(nChannels) ':'], 'CapSize', 10, 'LineWidth', nChannels_widths(nChannels)); hold on;
end
set(gca, 'YTick', [0.4 0.6 0.8], 'FontSize', 12);
set(gca, 'XTick', taus, 'FontSize', 12);
y = ylabel('r', 'FontSize', 15, 'rotation', 90);
xlabel('\tau', 'FontSize', 15);
axis([2 18 0.37 0.81]);
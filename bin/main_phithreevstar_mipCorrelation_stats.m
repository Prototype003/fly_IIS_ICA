%% Description

%{
Stats and figures for correlations between the phi3 and phistar values
across partition schemes
%}

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

%% Plot correlations at one parameter combination

nChannels = 2;
fly = 13;
tau = 3;
condition = 1;

figure;
bar(correlations{nChannels}.trial_averaged.phi_v_mistar(:, fly, condition, tau));
b = correlations{nChannels}.trial_averaged.phi_v_mistar(:, fly, condition, tau);

%% Portion of significant correlations within fly

portions = zeros(length(set_sizes), length(flies), length(conditions), length(taus));
portion_metric = 'phi_v_mistar';
q = 0.05;

for nChannels = 1 : length(correlations)
    for fly = 1 : length(flies)
        for matrix = 1 : length(correlation_matrices)
            for condition = 1 : length(conditions)
                for tau = 1 : length(taus)
                    
                    % Pre-corrected p values
                    p_mat = correlations{nChannels}.trial_averaged.([portion_metric p_suffix])(:, fly, condition, tau);
                    
                    % Get FDR threshold
                    [threshold_para, threshold_nonpara] = FDR(p_mat, q);
                    
                    % Find portion of p values which are a below the
                    % threshold
                    if numel(threshold_para) == 0 % Nothing is significant
                        portions(nChannels, fly, condition, tau) = 0;
                    else
                        % Portion of correlations which meet significance
                        sig_mat = p_mat < threshold_para;
                        
                        % Get and store portion of significant correlations
                        portions(nChannels, fly, condition, tau) = sum(sig_mat(:)) / numel(sig_mat);
                    end
                    
                end
            end
        end
    end
end

%% Replace non-sig correlations with r = 0

sig_suffix = '_sigs';

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
nChannels_widths = [1 2];
nChannels_markers = 'ox';
nChannels_offsets = [-0.1 0.1];

figure;
for nChannels = 1 : length(set_sizes)
    %subplot(1, length(nChannels), nChannels);
    %bar(correlations{nChannels}.trial_averaged.([metric mean_suffix])); hold on;
    
    errorbar(taus+nChannels_offsets(nChannels), correlations{nChannels}.trial_averaged.([metric mean_suffix]), correlations{nChannels}.trial_averaged.([metric std_suffix]), ['k' nChannels_markers(nChannels) '--'], 'CapSize', 10, 'LineWidth', nChannels_widths(nChannels)); hold on;
    
    %title([num2str(nChannels+2) ' channels']);
    ylabel('r');
    xlabel('tau');
    %set(gca, 'XTickLabel', taus);
end
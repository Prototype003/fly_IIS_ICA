%% Description

%{
Stats for correlations between the phi3 and phistar values
across partition schemes
%}

%% Setup

correlation_type = 'phi_v_phistar';

set_sizes = (3:4);
flies = (1:13);
conditions = (1);
taus = [4 8 16];

conditionSplit = 1;

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

% sig_suffix = '_sigs';
% 
% % Get logical indexes of significant correlations
% for nChannels = 1 : length(correlations)
%     for metric_counter = 1 : length(correlation_matrices)
%         metric = correlation_matrices{metric_counter};
%         
%         correlations{nChannels}.([metric sig_suffix]) = zeros(size(correlations{nChannels}.(metric)));
%         correlations{nChannels}.trial_averaged.([metric sig_suffix]) = zeros(size(correlations{nChannels}.(metric)));
%         
%         % FDR within each fly and tau
%         for fly = 1 : length(flies)
%             for condition = 1 : size(correlations{nChannels}.(metric), 3)
%                 for tau = 1 : length(taus)
%                     p_vector = correlations{nChannels}.([metric p_suffix])(:, fly, condition, tau);
%                     [threshold_para, threshold_nonpara] = FDR(p_vector, q);
%                     if numel(threshold_para) == 0
%                         threshold_para = 0;
%                     end
%                     correlations{nChannels}.([metric sig_suffix])(:, fly, condition, tau) = correlations{nChannels}.([metric p_suffix])(:, fly, condition, tau) < threshold_para;
%                     
%                     % repeat for trial averaged
%                     p_vector = correlations{nChannels}.trial_averaged.([metric p_suffix])(:, fly, condition, tau);
%                     [threshold_para, threshold_nonpara] = FDR(p_vector, q);
%                     if numel(threshold_para) == 0
%                         threshold_para = 0;
%                     end
%                     correlations{nChannels}.trial_averaged.([metric sig_suffix])(:, fly, condition, tau) = correlations{nChannels}.trial_averaged.([metric p_suffix])(:, fly, condition, tau) < threshold_para;
%                 end
%             end
%         end
%         
%     end
% end
% 
% % Replace non-significant correlations with 0
% for nChannels = 1 : length(correlations)
%     for metric_counter = 1 : length(correlation_matrices)
%         metric = correlation_matrices{metric_counter};
%         
%         for element = 1 : numel(correlations{nChannels}.(metric))
%             if correlations{nChannels}.([metric sig_suffix])(element) == 0
%                 correlations{nChannels}.(metric)(element) = 0;
%             end
%         end
%         for element = 1 : numel(correlations{nChannels}.trial_averaged.(metric))
%             if correlations{nChannels}.trial_averaged.([metric sig_suffix])(element) == 0
%                 correlations{nChannels}.trial_averaged.(metric)(element) = 0;
%             end
%         end
%         
%     end
% end

%% Build table for LME

disp('Building table');

% Determine total number of parameter combinations
table_length = 0;
for nChannels_counter = 1 : length(set_sizes)
    parameter_combos = numel(correlations{nChannels_counter}.(correlation_type));
    table_length = table_length + parameter_combos;
end

% Extract all values and parameters into vectors
table_raw = struct();

table_raw.nChannels = zeros(table_length, 1);
table_raw.condition = zeros(table_length, 1);
table_raw.fly = zeros(table_length, 1);
table_raw.tau = zeros(table_length, 1);
table_raw.set = zeros(table_length, 1);
table_raw.correlation = zeros(table_length, 1);

row_counter = 1;
set_label_starts = [1 106 456];
for nChannels_counter = 1 : length(correlations)
    for condition = 1 : size(correlations{nChannels_counter}.(correlation_type), 3)
        for fly = 1 : size(correlations{nChannels_counter}.(correlation_type), 2)
            for tau_counter = 1 : size(correlations{nChannels_counter}.(correlation_type), 4)
                set_label = set_label_starts(nChannels_counter);
                for set = 1 : size(correlations{nChannels_counter}.(correlation_type), 1)
                    table_raw.nChannels(row_counter) = correlations{nChannels_counter}.nChannels;
                    table_raw.condition(row_counter) = condition;
                    table_raw.fly(row_counter) = fly;
                    table_raw.tau(row_counter) = correlations{nChannels_counter}.taus(tau_counter);
                    table_raw.set(row_counter) = set_label; set_label = set_label + 1;
                    table_raw.correlation(row_counter) = fisher_rz(correlations{nChannels_counter}.trial_averaged.(correlation_type)(set, fly, 1, tau_counter));
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end

table_headings = {'fly', 'condition', 'nChannels', 'tau', 'set', 'correlation'};
correlation_table = table(table_raw.fly, table_raw.condition, table_raw.nChannels, table_raw.tau, table_raw.set, table_raw.correlation, 'VariableNames', table_headings);

%% LME Model

% We want to nest set within fly

model_spec = 'correlation ~ condition + nChannels + tau + (1|fly) + (1|fly:set)';

disp(['fitting model: ' model_spec]);

model_full = fitlme(correlation_table, model_spec);

disp('full model built');

%% Null Models

model_null_specs = {...
    'correlation ~ nChannels + tau + (1|fly) + (1|fly:set)',... % condition null model
    'correlation ~ condition + tau + (1|fly) + (1|fly:set)',... % nChannels null model
    'correlation ~ condition + nChannels + (1|fly) + (1|fly:set)'... % tau null model
    };

% Build null models
model_nulls = cell(length(model_null_specs), 1);
for null = 1 : length(model_null_specs)
    disp(['fitting null model: ' model_null_specs{null}]);
    model_nulls{null} = fitlme(correlation_table, model_null_specs{null});
    disp('null model built');
end

%% Likelihood ratio tests for main effects

compare(model_nulls{1}, model_full);
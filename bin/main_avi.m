%% DESCRIPTION

%{
This script tests the first hypothesis: phi in the air condition is greater than phi in the isoflurane condition

Tests for all nChannels in loaded datafile

Test is conducted using LME with model comparisons (to the null)

%}

%% SETUP

phi_type = 'three';

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_phi' phi_type...
    ];

results_directory = 'analysis_results/';
results_filename = [data_filename '_lmeStats'];

%% LOAD

disp('loading');

load([data_directory data_filename '.mat']);

disp('loaded');

%% Preprocess

for nChannels_counter = 1 : length(phis)
    % So I don't need to keep on typing ...}.phi_threes
    phis{nChannels_counter}.phis = phis{nChannels_counter}.(['phi_' phi_type 's']);
    
    % Average across trials
    phis{nChannels_counter}.phis = mean(phis{nChannels_counter}.phis, 2);
end

%% Build table for LME

disp('building table');

% Determine total number of parameter combinations
table_length = 0;
for nChannels_counter = 1 : length(phis)
    parameter_combos = numel(phis{nChannels_counter}.phis);
    table_length = table_length + parameter_combos;
end

% Extract all values and parameters into vectors
table_raw = struct();
table_raw.nChannels = zeros(table_length, 1);
table_raw.set = zeros(table_length, 1);
table_raw.trial = zeros(table_length, 1);
table_raw.fly = zeros(table_length, 1);
table_raw.condition = zeros(table_length, 1);
table_raw.tau = zeros(table_length, 1);
table_raw.phi = zeros(table_length, 1);
row_counter = 1;
% Loop is based on natural ordering: fly, condition, tau, nChannels, set, trial
% Assumes same number of flies, conditions and taus for all nChannels
for fly = 1 : size(phis{1}.phis, 3)
    for condition = 1 : size(phis{1}.phis, 4)
        for tau = 1 : size(phis{1}.phis, 5)
            for nChannels = 1 : length(phis)
                for set = 1 : size(phis{nChannels}.phis, 1)
                    for trial = 1 : size(phis{nChannels}.phis, 2)
                        table_raw.nChannels(row_counter) = phis{nChannels}.nChannels;
                        table_raw.set(row_counter) = set;
                        table_raw.trial(row_counter) = trial;
                        table_raw.fly(row_counter) = fly;
                        table_raw.condition(row_counter) = condition;
                        table_raw.tau(row_counter) = tau;
                        table_raw.phi(row_counter) = log(phis{nChannels}.phis(set, trial, fly, condition, tau));
                        row_counter = row_counter + 1;
                    end
                end
            end
        end
    end
end

% Natural ordering: fly, condition, tau, nChannels, set, trial
table_headings = {'fly', 'condition', 'tau', 'nChannels', 'set', 'trial', 'phi'};
phi_table = table(table_raw.fly, table_raw.condition, table_raw.tau, table_raw.nChannels, table_raw.set, table_raw.trial, table_raw.phi, 'VariableNames', table_headings);

disp('table built');

%% LME model

% We want to take everything into account
% channel sets are nested within fly
% channel sets are crossed across trials
% Main effects: set size (nChannels), condition, tau
% Nested effects: (1|fly/set)=(1|fly) +  OR (1|fly/trial/set)
% Crossed effects: (1|trial) (may just ignore this)
% Consequences of not include trial in the model? Less fit, so test is more conservative?

model_spec = 'phi ~ nChannels + condition + tau + (1|fly) + (1|fly:set)';


disp(['fitting model: ' model_spec]);

model_full = fitlme(phi_table, model_spec);

disp('full model built');

%% Build null models for comparison

model_null_specs = {...
    'phi ~ condition + tau + (1|fly) + (1|fly:set)',... % nChannels null model
    'phi ~ nChannels + tau + (1|fly) + (1|fly:set)',... % condition null model
    'phi ~ nChannels + condition + (1|fly) + (1|fly:set)'... % tau null model
    };

% Build null models
model_nulls = cell(length(model_null_specs), 1);
for null = 1 : length(model_null_specs)
    disp(['fitting null model: ' model_null_specs{null}]);
    model_nulls{null} = fitlme(phi_table, model_null_specs{null});
    disp('null model built');
end

%% Likelihood ratio test
% % Likelihood ratio test is anticonservative when testing for fixed effects, so use simulated test
% 
% iterations = 1000;
% 
% options = statset('LinearMixedModel');
% options.UseParallel = true;
% 
% compute_pool = parpool();
% 
% model_comparisons = cell(length(model_nulls), 1);
% model_comparisons_sim = cell(length(model_nulls), 1);
% for null = 1 : length(model_nulls)
%     disp(['conducting likelihood ratio tests (sim) on: ' model_null_specs{null}]);
%     [model_comparisons{null}, model_comparisons_sim{null}] = compare(model_nulls{null}, model_full, 'NSim', iterations, 'Options', options);
%     %model_comparisons{null} = compare(model_nulls{null}, model_full);
%     disp('completed');
% end
% 
% delete(compute_pool);

%% Save results

% disp('Saving');
% if ~isdir(results_directory)
%     mkdir(results_directory)
% end
% save([results_directory results_filename '.mat'],...
%     'phi_table',...
%     'model_spec',...
%     'model_full',...
%     'model_null_specs',...
%     'model_nulls',...
%     'model_comparisons',...
%     'model_comparisons_sim'...
%     );
% 
% disp('Saved');
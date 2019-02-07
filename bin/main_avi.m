%% DESCRIPTION

%{
This script tests the first hypothesis: phi in the air condition is greater than phi in the isoflurane condition

Tests for all nChannels in loaded datafile

Test is conducted using LME with model comparisons (to the null)

%}

%% SETUP

bin_location = './';
phi_type = 'phi_three'; % 'phi_three' or 'phi_star_gaussian'
global_tpm = 0;

results_directory = 'analysis_results/';
results_filename = [phi_type '_globalTPM' num2str(global_tpm) '_lmeStats'];

%% LOAD

addpath('figure_code/');
[phis, measure_strings{1}] = phi_load(phi_type, global_tpm, bin_location);

%% Preprocess

for nChannels_counter = 1 : length(phis)
    % Average across trials
    phis{nChannels_counter}.phis_orig = phis{nChannels_counter}.phis; % Backup phi values
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
table_raw.center = zeros(table_length, 1);
table_raw.distance = zeros(table_length, 1);
table_raw.phi = zeros(table_length, 1);
row_counter = 1;
set_label_starts = [1 106 456]; % E.g. set 1 of 2 channels is given a different label to set 1 of 4 channels
tau_labels = [4 8 16];
% Loop is based on natural ordering: fly, condition, tau, nChannels, set, trial
% Assumes same number of flies, conditions and taus for all nChannels
for fly = 1 : size(phis{1}.phis, 3)
    for condition = 1 : size(phis{1}.phis, 4)
        for tau = 1 : size(phis{1}.phis, 5)
            for nChannels = 1 : length(phis)
                channel_sets = double(phis{nChannels}.channel_sets);
                centers = mean(channel_sets, 2); % Mean across channels in each set
                distances = channel_set_distances(channel_sets);
                set_label = set_label_starts(nChannels_counter);
                for set_counter = 1 : size(phis{nChannels}.phis, 1)
                    for trial = 1 : size(phis{nChannels}.phis, 2)
                        table_raw.nChannels(row_counter) = phis{nChannels}.nChannels;
                        table_raw.set(row_counter) = set_counter;
                        table_raw.trial(row_counter) = trial;
                        table_raw.fly(row_counter) = fly;
                        table_raw.condition(row_counter) = condition;
                        table_raw.tau(row_counter) = phis{nChannels}.taus(tau);
%                         if phis{nChannels}.phis(set_counter, trial, fly, condition, tau) <= 0
%                             tmp = phis{nChannels}.phis(:);
%                             tmp(tmp <= 0) = inf;
%                             table_raw.phi(row_counter) = log(min(tmp));
%                         else
%                             table_raw.phi(row_counter) = log(phis{nChannels}.phis(set_counter, trial, fly, condition, tau));
%                         end
                        table_raw.center(row_counter) = centers(set_counter);
                        table_raw.distance(row_counter) = distances(set_counter);
                        table_raw.phi(row_counter) = phis{nChannels}.phis(set_counter, trial, fly, condition, tau);
                        row_counter = row_counter + 1;
                    end
                    set_label = set_label + 1;
                end
            end
        end
    end
end

% Apply dummy coding to condition, set-size, tau

% Natural ordering: fly, condition, tau, nChannels, set, trial
table_headings = {'fly', 'condition', 'tau', 'nChannels', 'set', 'center', 'distance', 'trial', 'phi'};
phi_table = table(table_raw.fly, table_raw.condition, table_raw.tau, (table_raw.nChannels), table_raw.set, table_raw.center, table_raw.distance, table_raw.trial, log(table_raw.phi), 'VariableNames', table_headings);
%phi_table.fly = nominal(phi_table.fly);
%phi_table.condition = nominal(phi_table.condition);
%phi_table.set = nominal(phi_table.set);
phi_dataset = table2dataset(phi_table);
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
%model_spec = 'phi ~ nChannels + tau + center*distance + (1|fly) + (1|fly:set)';
model_spec = 'phi ~ nChannels + condition + tau + (1|fly)';

model_spec = 'phi ~ condition + nChannels + tau + (1+condition+nChannels+tau|fly)';

disp(['fitting model: ' model_spec]);

% Warning - GLME using Laplace method (needed for likelihood ratio test) takes a long time
% Also, GLME using lognormal function still has heteroscedasticity

model_full = fitlme(phi_table, model_spec);
%model_full = fitlme(phi_table(phi_table.condition==1, :), model_spec);
%model_full = fitglme(phi_table(phi_table.fly==1, :), model_spec, 'Distribution', 'Normal', 'Link', 'log', 'FitMethod', 'Laplace'); % Log-normal distribution has normal distribution and log link function

disp('full model built');

figure; scatter(phi_table.phi, fitted(model_full));
%figure; scatter(phi_table.phi(phi_table.fly==1, :), fitted(model_full)); line([0 0.35], [0 0.35]);

%% Build null models for comparison

model_null_specs = {...
    'phi ~ nChannels + tau + (1|fly) + (1|fly:set)',... % condition null model
    'phi ~ condition + tau + (1|fly) + (1|fly:set)',... % nChannels null model
    'phi ~ nChannels + condition + (1|fly) + (1|fly:set)'... % tau null model
    'phi ~ nChannels + condition + tau + (1|fly) + (1|fly:set)'... % full main effect null model
    };

model_null_specs = {...
    'phi ~ nChannels + tau + (1+condition+nChannels+tau|fly)',... % condition null model
    'phi ~ condition + tau + (1+condition+nChannels+tau|fly)',... % nChannels null model
    'phi ~ nChannels + condition + (1+condition+nChannels+tau|fly)'... % tau null model
    'phi ~ nChannels + condition + tau + (1+condition+nChannels+tau|fly)'... % full main effect null model
    };

% model_null_specs = {...
%     'phi ~ nChannels + tau + center*distance + (1|fly) + (1|fly:set)',... % condition null model
%     'phi ~ tau + center*distance + (1|fly) + (1|fly:set)',... % nChannels null model
%     'phi ~ nChannels + center*distance + (1|fly) + (1|fly:set)'... % tau null model
%     'phi ~ nChannels + tau + center + distance + (1|fly) + (1|fly:set)'... % center:distance interaction null model
%     };

model_null_specs = {'phi ~ nChannels + tau + (1+condition+nChannels+tau|fly)'};

% Build null models
model_nulls = cell(length(model_null_specs), 1);
for null = 1 : length(model_null_specs)
    disp(['fitting null model: ' model_null_specs{null}]);
    model_nulls{null} = fitlme(phi_table, model_null_specs{null});
    disp('null model built');
end

%phis{nChannels_counter}.phis = phis{nChannels_counter}.phis_orig; % restore backed up (non-averaged) values

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

% Save results

disp('Saving');
if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename '.mat'],...
    'phi_table',...
    'model_spec',...
    'model_full',...
    'model_null_specs',...
    'model_nulls'...
    );
    %'model_comparisons',...
    %'model_comparisons_sim'...
    %);

disp('Saved');

%% Post-hoc tests : condition

% Average across channel sets to get single value per fly
sigs = zeros(length(phis), length(phis{1}.taus));
zs = zeros(length(phis), length(phis{1}.taus));
reduced = zeros(length(phis), size(phis{1}.phis, 3), size(phis{1}.phis, 5));
for nChannels = 1 : length(phis)
    phis{nChannels}.phis_channelAveraged = permute(mean(phis{nChannels}.phis, 1), [3 4 5 1 2]);
    reduced(nChannels, :, :) = permute(...
        phis{nChannels}.phis_channelAveraged(:, 1, :) > phis{nChannels}.phis_channelAveraged(:, 2, :),...
        [1 3 2]...
        );
    
    % Test if reduction in phi is significant
    for tau = 1 : size(phis{nChannels}.phis_channelAveraged, 3)
        air = phis{nChannels}.phis_channelAveraged(:, 1, tau);
        iso = phis{nChannels}.phis_channelAveraged(:, 2, tau);
        [sigs(nChannels, tau), decision, stats] = signrank(air, iso, 'Tail', 'right', 'Method', 'approx');
        zs(nChannels, tau) = stats.zval;
    end
end

% Find which channel sets have a significant decrease in phi across flies (using trial-averaged values)
q = 0.05;
sig_portions = zeros(length(phis), length(phis{1}.taus));
sig_counts = zeros(length(phis), length(phis{1}.taus));
for nChannels = 1 : length(phis)
    phis{nChannels}.channel_sigs = zeros(size(phis{nChannels}.phis, 1), size(phis{nChannels}.phis, 5));
    phis{nChannels}.channel_sigs_fdr = zeros(size(phis{nChannels}.phis, 1), size(phis{nChannels}.phis, 5));
    for tau = 1 : size(phis{nChannels}.phis, 5)
        
        for channel_set = 1 : size(phis{nChannels}.phis, 1)
            
            air = squeeze(phis{nChannels}.phis(channel_set, 1, :, 1, tau));
            iso = squeeze(phis{nChannels}.phis(channel_set, 1, :, 2, tau));
            
            % Conduct test across flies
            [phis{nChannels}.channel_sigs(channel_set, tau), decision] =...
                signrank(air(:), iso(:), 'Tail', 'right'); % paired non-parametric test
            
        end
        
        % FDR correction
        phis{nChannels}.channel_sigs_fdr(:, tau) = fdr_correct(phis{nChannels}.channel_sigs(:, tau), q);
        
    end
    
    sig_portions(nChannels, :) = mean(phis{nChannels}.channel_sigs_fdr, 1);
    sig_counts(nChannels, :) = sum(phis{nChannels}.channel_sigs_fdr, 1);
end

% 

%% Post-hoc tests : nChannels

% Average across sets, conditions, and taus
values = zeros(length(phis), size(phis{1}.phis, 3));
for nChannels = 1 : length(phis)
    values(nChannels, :) = permute(mean(mean(mean(phis{nChannels}.phis, 1), 4), 5), [3 1 2 4 5]);
end

% Conduct tests
clear sigs
test_counter = 1;
for nChannels_1 = 1 : length(phis)
    for nChannels_2 = 2 : length(phis)
        if nChannels_1 ~= nChannels_2 && nChannels_2 > nChannels_1
            [p, decision, stats] = signrank(values(nChannels_1, :), values(nChannels_2, :), 'Tail', 'left', 'Method', 'approx');
            disp([num2str(nChannels_1) ' ' num2str(nChannels_2) ': Z = ' num2str(stats.zval) ' p = ' num2str(p)]);
            sigs(test_counter) = p;
            test_counter = test_counter + 1;
        end
    end
end

%% Post-hoc tests : tau

% Average across all sets (from all set sizes), and conditions
values = zeros(length(phis{1}.taus), size(phis{1}.phis, 3));
set_counter = 0;
for nChannels = 1 : length(phis)
    set_counter = set_counter + size(phis{nChannels}.phis, 1)*size(phis{nChannels}.phis, 4);
    values(:, :) = values(:, :) + permute(sum(sum(phis{nChannels}.phis, 1), 4), [5 3 1 2 4]);
end
values = values ./ set_counter;

% Conduct tests
clear sigs
test_counter = 1;
for nChannels_1 = 1 : length(phis)
    for nChannels_2 = 2 : length(phis)
        if nChannels_1 ~= nChannels_2 && nChannels_2 > nChannels_1
            [p, decision, stats] = signrank(values(nChannels_1, :), values(nChannels_2, :), 'Method', 'approx');
            disp([num2str(nChannels_1) ' ' num2str(nChannels_2) ': Z = ' num2str(stats.zval) ' p = ' num2str(p)]);
            sigs(test_counter) = p;
            test_counter = test_counter + 1;
        end
    end
end
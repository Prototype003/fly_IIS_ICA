%% DESCRIPTION

%{
Correlates phi3 and phistar (I-star) values at each partition

Reasoning behind correlating phi3 and I-star values:
    MIP = partition which has smallest difference in information to full
        system
    Phi-3 = partition with smallest phi
    Phistar MIP = partition with smallest normalised phi
    Phistar = I - I* ; so MIP should have largest I*
    Phi-3 aims to minimise phi, phistar aims to maximise I*, so negative
    correlation between phi3 and I*

Correlate at each channel set, across states?

%}

%% SETUP
flies = (1:13);
nChannels_set = (3:4);
possible_bipartitions = [3, 7];
conditions = (1:2);
taus = [4 8 16];
phithree_fly_separation = 4; % Results are split by fly for nChannels == 4

conditionSplit = 1;

data_directory = 'results/preformatted_results/';
phistar_file_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_medianSplit0_phistar_allPartitionsfly';
%phithree_file_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_allPartitionsfly5';
suffix = '.mat';

results_directory = 'analysis_results/';
results_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels' num2str(nChannels_set(1)) 't' num2str(nChannels_set(2)) '_phithreevstar_correlations_conditionSplit' num2str(conditionSplit) suffix];

%% Correlate per fly

%{
Data structure:
Correlations at each channel set, at each trial
Average correlation coefficients across trials
Average correltion coefficients across channels

Dimension reminder:
state_partitions_phis = (states x partitions x sets x flies x conditions x taus)
phi-star = (partitions x sets x trials x flies x conditions x taus)
%}

%%

if conditionSplit == 1
    
    % Conduct analysis at the per fly level
    fly_correlations = cell(1, length(flies));
    for fly_counter = 1 : length(flies)
        fly = flies(fly_counter);
        disp(['fly ' num2str(fly)]);
        
        % Load phi-3 results
        disp('Loading phi-3');
        phi_threes = load_nChannels(fly, nChannels_set, phithree_fly_separation);
        disp('phi-3 loaded');
        
        % Load phi-* results
        disp('Loading phi-*');
        phi_stars = load([data_directory phistar_file_prefix num2str(fly) suffix]);
        phi_stars = phi_stars.phis;
        phi_stars = phi_stars(nChannels_set-1);
        disp('phi-* loaded');
        
        % Conduct analysis at the nChannels level
        % Results storage structure
        correlations = cell(1, length(nChannels_set));
        for nChannels_counter = 1 : length(nChannels_set)
            nChannels = nChannels_set(nChannels_counter);
            
            nPartitions = possible_bipartitions(nChannels_counter);
            
            % We want to calculate correlation for each channel set
            nSets = size(phi_threes{nChannels_counter}.phi_threes, 1);
            
            % Correlation storage matrices
            phi_v_phistar = zeros(nSets, 1, length(conditions), length(taus));
            phi_v_phistar_ps = zeros(nSets, 1, length(conditions), length(taus));
            phi_v_mi = zeros(nSets, 1, length(conditions), length(taus));
            phi_v_mi_ps = zeros(nSets, 1, length(conditions), length(taus));
            phi_v_mistar = zeros(nSets, 1, length(conditions), length(taus));
            phi_v_mistar_ps = zeros(nSets, 1, length(conditions), length(taus));
            
            trial_averaged = struct();
            trial_averaged.phi_v_phistar = zeros(nSets, 1, length(conditions), length(taus));
            trial_averaged.phi_v_phistar_ps = zeros(nSets, 1, length(conditions), length(taus));
            trial_averaged.phi_v_mi = zeros(nSets, 1, length(conditions), length(taus));
            trial_averaged.phi_v_mi_ps = zeros(nSets, 1, length(conditions), length(taus));
            trial_averaged.phi_v_mistar = zeros(nSets, 1, length(conditions), length(taus));
            trial_averaged.phi_v_mistar_ps = zeros(nSets, 1, length(conditions), length(taus));
            
            for channel_set = 1 : nSets
                
                for condition = 1 : length(conditions)
                    for tau_counter = 1 : length(taus)
                        
                        % First get trial-averaged phi-star values (this will be same
                        % for all phi-3 states), and get partition order
                        phi_star_partition_phistars = zeros(size(phi_stars{nChannels_counter}.partitions, 3), 2*nPartitions); % trials x bidirectional partitions
                        phi_star_partition_mis = zeros(size(phi_stars{nChannels_counter}.partitions, 3), 2*nPartitions);
                        phi_star_partition_mistars = zeros(size(phi_stars{nChannels_counter}.partitions, 3), 2*nPartitions);
                        value_counter = 1; % Used to keep track of phi_star_partition_xxx position, as partition_counter increase by 1, value_counter increases by 2
                        partition_order = cell(1, nPartitions);
                        partition_position = 1;
                        for partition_counter = 1 : size(phi_stars{nChannels_counter}.partitions, 1)
                            partition = phi_stars{nChannels_counter}.partitions{partition_counter, channel_set, 1, 1, condition, tau_counter};
                            if bipartition_check(partition)
                                % Note that phi-star calculation was not
                                % parallelised, thus partition order is
                                % guaranteed to be the same across trials and
                                % other parameters
                                
                                % Add partition to partition order
                                partition_order{partition_position} = partition;
                                partition_position = partition_position + 1;
                                
                                % Get values across trials
                                phistars = squeeze(phi_stars{nChannels_counter}.partitions_phi_stars(partition_counter, channel_set, :, 1, condition, tau_counter));
                                mis = squeeze(phi_stars{nChannels_counter}.partitions_mis(partition_counter, channel_set, :, 1, condition, tau_counter));
                                mi_stars = squeeze(phi_stars{nChannels_counter}.partitions_mi_stars(partition_counter, channel_set, :, 1, condition, tau_counter));
                                
                                % Add trial values (x2)
                                phi_star_partition_phistars(:, value_counter) = phistars;
                                phi_star_partition_mis(:, value_counter) = mis;
                                phi_star_partition_mistars(:, value_counter) = mi_stars;
                                value_counter = value_counter + 1;
                                phi_star_partition_phistars(:, value_counter) = phistars;
                                phi_star_partition_mis(:, value_counter) = mis;
                                phi_star_partition_mistars(:, value_counter) = mi_stars;
                                value_counter = value_counter + 1;
                                
                            end
                        end
                        
                        % Get phi3 values for each partition, in the same
                        % order as for phistar
                        phi_three_partition_phis = zeros(size(phi_threes{nChannels_counter}.state_counters, 3), 2*nPartitions); % trials x directional partitions
                        value_counter = 1;
                        for partition_counter = 1 : length(partition_order)
                            partition = partition_order{partition_counter};
                            
                            state_phis = zeros(2, size(phi_threes{nChannels_counter}.state_partitions, 1));
                            
                            % Get values for all states (average across states)
                            for state = 1 : size(phi_threes{nChannels_counter}.state_partitions, 1)
                                
                                % Find partitions which are the same
                                partition_matches = 0;
                                for phi_three_partition_counter = 1 : size(phi_threes{nChannels_counter}.state_partitions, 2)
                                    phi_three_partition = phi_threes{nChannels_counter}.state_partitions{state, phi_three_partition_counter, channel_set, 1, condition, tau_counter};
                                    
                                    if mip_equal(partition, relabel(phi_three_partition))
                                        % Record phi value for the partition
                                        state_phis(partition_matches+1, state) = phi_threes{nChannels_counter}.state_partitions_phis(state, phi_three_partition_counter, channel_set, 1, condition, tau_counter);
                                        partition_matches = partition_matches + 1;
                                    end
                                    
                                    if partition_matches >= 2
                                        break % There are only 2 phi3 partitions for each phi* partition (due to bidirectionality)
                                    end
                                end
                                
                            end
                            
                            % Get weighted phi average at each trial
                            trial_phis = zeros(size(phi_threes{nChannels_counter}.state_counters, 3), 2);
                            for trial = 1 : size(phi_threes{nChannels_counter}.state_counters, 3)
                                % Get counts for each state in the trial
                                state_counts = phi_threes{nChannels_counter}.state_counters(:, channel_set, trial, 1, condition, tau_counter);
                                
                                % Get weighted average for the trial
                                trial_state_phis = zeros(size(state_phis));
                                trial_state_phis(1, :) = state_phis(1, :) .* state_counts';
                                trial_state_phis(2, :) = state_phis(2, :) .* state_counts';
                                trial_state_phis = sum(trial_state_phis, 2) / sum(state_counts);
                                
                                trial_phis(trial, :) = trial_state_phis;
                            end
                            
                            % Store trial values
                            phi_three_partition_phis(:, value_counter:value_counter+1) = trial_phis;
                            value_counter = value_counter + 2;
                            
                        end
                        
                        % Correlate phi3 and phi* values at each state, using
                        % all trials
                        [r, p] = corrcoef(phi_star_partition_phistars(:), phi_three_partition_phis(:));
                        phi_v_phistar(channel_set, 1, condition, tau_counter) = r(1, 2);
                        phi_v_phistar_ps(channel_set, 1, condition, tau_counter) = p(1, 2);
                        
                        [r, p] = corrcoef(phi_star_partition_mis(:), phi_three_partition_phis(:));
                        phi_v_mi(channel_set, 1, condition, tau_counter) = r(1, 2);
                        phi_v_mi_ps(channel_set, 1, condition, tau_counter) = p(1, 2);
                        
                        [r, p] = corrcoef(phi_star_partition_mistars(:), phi_three_partition_phis(:));
                        phi_v_mistar(channel_set, 1, condition, tau_counter) = r(1, 2);
                        phi_v_mistar_ps(channel_set, 1, condition, tau_counter) = p(1, 2);
                        
                        % correlate after average across trials
                        [r, p] = corrcoef(mean(phi_star_partition_phistars, 1), mean(phi_three_partition_phis, 1));
                        trial_averaged.phi_v_phistar(channel_set, 1, condition, tau_counter) = r(1, 2);
                        trial_averaged.phi_v_phistar_ps(channel_set, 1, condition, tau_counter) = p(1, 2);
                        
                        [r, p] = corrcoef(mean(phi_star_partition_mis, 1), mean(phi_three_partition_phis, 1));
                        trial_averaged.phi_v_mi(channel_set, 1, condition, tau_counter) = r(1, 2);
                        trial_averaged.phi_v_mi_ps(channel_set, 1, condition, tau_counter) = p(1, 2);
                        
                        [r, p] = corrcoef(mean(phi_star_partition_mistars, 1), mean(phi_three_partition_phis, 1));
                        trial_averaged.phi_v_mistar(channel_set, 1, condition, tau_counter) = r(1, 2);
                        trial_averaged.phi_v_mistar_ps(channel_set, 1, condition, tau_counter) = p(1, 2);
                    end
                end
                
            end
            
            % Store
            correlations{nChannels_counter} = struct();
            correlations{nChannels_counter}.phi_v_phistar = phi_v_phistar;
            correlations{nChannels_counter}.phi_v_phistar_ps = phi_v_phistar_ps;
            correlations{nChannels_counter}.phi_v_mi = phi_v_mi;
            correlations{nChannels_counter}.phi_v_mi_ps = phi_v_mi_ps;
            correlations{nChannels_counter}.phi_v_mistar = phi_v_mistar;
            correlations{nChannels_counter}.phi_v_mistar_ps = phi_v_mistar_ps;
            correlations{nChannels_counter}.trial_averaged = trial_averaged;
            
        end
        
        % Store
        fly_correlations{fly_counter} = correlations;
    end
    
end

%%

if conditionSplit == 0
    
    % Conduct analysis at the per fly level
    fly_correlations = cell(1, length(flies));
    for fly_counter = 1 : length(flies)
        fly = flies(fly_counter);
        disp(['fly ' num2str(fly)]);
        
        % Load phi-3 results
        disp('Loading phi-3');
        phi_threes = load_nChannels(fly, nChannels_set, phithree_fly_separation);
        disp('phi-3 loaded');
        
        % Load phi-* results
        disp('Loading phi-*');
        phi_stars = load([data_directory phistar_file_prefix num2str(fly) suffix]);
        phi_stars = phi_stars.phis;
        phi_stars = phi_stars(nChannels_set-1);
        disp('phi-* loaded');
        
        % Conduct analysis at the nChannels level
        % Results storage structure
        correlations = cell(1, length(nChannels_set));
        for nChannels_counter = 1 : length(nChannels_set)
            nChannels = nChannels_set(nChannels_counter);
            
            nPartitions = possible_bipartitions(nChannels_counter);
            
            % We want to calculate correlation for each channel set
            nSets = size(phi_threes{nChannels_counter}.phi_threes, 1);
            
            % Correlation storage matrices
            phi_v_phistar = zeros(nSets, 1, length(taus));
            phi_v_phistar_ps = zeros(nSets, 1, length(taus));
            phi_v_mi = zeros(nSets, 1, length(taus));
            phi_v_mi_ps = zeros(nSets, 1, length(taus));
            phi_v_mistar = zeros(nSets, 1, length(taus));
            phi_v_mistar_ps = zeros(nSets, 1, length(taus));
            
            trial_averaged = struct();
            trial_averaged.phi_v_phistar = zeros(nSets, 1, length(taus));
            trial_averaged.phi_v_phistar_ps = zeros(nSets, 1, length(taus));
            trial_averaged.phi_v_mi = zeros(nSets, 1, length(taus));
            trial_averaged.phi_v_mi_ps = zeros(nSets, 1, length(taus));
            trial_averaged.phi_v_mistar = zeros(nSets, 1, length(taus));
            trial_averaged.phi_v_mistar_ps = zeros(nSets, 1, length(taus));
            
            for channel_set = 1 : nSets
                
                for tau_counter = 1 : length(taus)
                    
                    % First get trial-averaged phi-star values (this will be same
                    % for all phi-3 states), and get partition order
                    phi_star_partition_phistars = zeros(size(phi_stars{nChannels_counter}.partitions, 3)*length(conditions) , 2*nPartitions); % trials*conditions x bidirectional partitions
                    phi_star_partition_mis = zeros(size(phi_stars{nChannels_counter}.partitions, 3)*length(conditions), 2*nPartitions);
                    phi_star_partition_mistars = zeros(size(phi_stars{nChannels_counter}.partitions, 3)*length(conditions), 2*nPartitions);
                    value_counter = 1; % Used to keep track of phi_star_partition_xxx position, as partition_counter increase by 1, value_counter increases by 2
                    partition_order = cell(1, nPartitions);
                    partition_position = 1;
                    for partition_counter = 1 : size(phi_stars{nChannels_counter}.partitions, 1)
                        partition = phi_stars{nChannels_counter}.partitions{partition_counter, channel_set, 1, 1, 1, tau_counter};
                        if bipartition_check(partition)
                            % Note that phi-star calculation was not
                            % parallelised, thus partition order is
                            % guaranteed to be the same across trials and
                            % other parameters
                            
                            % Add partition to partition order
                            partition_order{partition_position} = partition;
                            partition_position = partition_position + 1;
                            
                            % Get values across trials
                            phistars = squeeze(phi_stars{nChannels_counter}.partitions_phi_stars(partition_counter, channel_set, :, 1, :, tau_counter));
                            mis = squeeze(phi_stars{nChannels_counter}.partitions_mis(partition_counter, channel_set, :, 1, :, tau_counter));
                            mi_stars = squeeze(phi_stars{nChannels_counter}.partitions_mi_stars(partition_counter, channel_set, :, 1, :, tau_counter));
                            
                            % Add trial values (x2)
                            phi_star_partition_phistars(:, value_counter) = phistars(:);
                            phi_star_partition_mis(:, value_counter) = mis(:);
                            phi_star_partition_mistars(:, value_counter) = mi_stars(:);
                            value_counter = value_counter + 1;
                            phi_star_partition_phistars(:, value_counter) = phistars(:);
                            phi_star_partition_mis(:, value_counter) = mis(:);
                            phi_star_partition_mistars(:, value_counter) = mi_stars(:);
                            value_counter = value_counter + 1;
                            
                        end
                    end
                    
                    % Get phi3 values for each partition, in the same
                    % order as for phistar
                    phi_three_partition_phis = zeros(size(phi_threes{nChannels_counter}.state_counters, 3)*length(conditions), 2*nPartitions); % trials*conditions x directional partitions
                    value_counter = 1;
                    for partition_counter = 1 : length(partition_order)
                        partition = partition_order{partition_counter};
                        
                        state_phis = zeros(2, size(phi_threes{nChannels_counter}.state_partitions, 1), length(conditions));
                        
                        % Get values for all states (average across states)
                        for state = 1 : size(phi_threes{nChannels_counter}.state_partitions, 1)
                            
                            for condition = 1 : length(conditions)
                                % Find partitions which are the same
                                partition_matches = 0;
                                for phi_three_partition_counter = 1 : size(phi_threes{nChannels_counter}.state_partitions, 2)
                                    phi_three_partition = phi_threes{nChannels_counter}.state_partitions{state, phi_three_partition_counter, channel_set, 1, condition, tau_counter};
                                    
                                    if mip_equal(partition, relabel(phi_three_partition))
                                        % Record phi value for the partition
                                        state_phis(partition_matches+1, state, condition) = phi_threes{nChannels_counter}.state_partitions_phis(state, phi_three_partition_counter, channel_set, 1, condition, tau_counter);
                                        partition_matches = partition_matches + 1;
                                    end
                                    
                                    if partition_matches >= 2
                                        break % There are only 2 phi3 partitions for each phi* partition (due to bidirectionality)
                                    end
                                end
                            end
                            
                        end
                        
                        % Get weighted phi average at each trial
                        trial_phis = zeros(size(phi_threes{nChannels_counter}.state_counters, 3)*length(conditions), 2); % trials*conditions x directional partitions
                        trial_condition = 1;
                        for condition = 1 : length(conditions)
                            for trial = 1 : size(phi_threes{nChannels_counter}.state_counters, 3) % trial in the inner loop means that trials will be grouped by condition, corresponding to phi* values
                                % Get counts for each state in the trial
                                state_counts = phi_threes{nChannels_counter}.state_counters(:, channel_set, trial, 1, condition, tau_counter);
                                
                                % Get weighted average for the trial
                                trial_state_phis = zeros(size(state_phis, 1), size(state_phis, 2));
                                trial_state_phis(1, :) = state_phis(1, :, condition) .* state_counts';
                                trial_state_phis(2, :) = state_phis(2, :, condition) .* state_counts';
                                trial_state_phis = sum(trial_state_phis, 2) / sum(state_counts);
                                
                                trial_phis(trial_condition, :) = trial_state_phis;
                                trial_condition = trial_condition + 1;
                            end
                        end
                        
                        % Store trial values
                        phi_three_partition_phis(:, value_counter:value_counter+1) = trial_phis;
                        value_counter = value_counter + 2;
                        
                    end
                    
                    % Correlate phi3 and phi* values at each state, using
                    % all trials
                    [r, p] = corrcoef(phi_star_partition_phistars(:), phi_three_partition_phis(:));
                    phi_v_phistar(channel_set, 1, tau_counter) = r(1, 2);
                    phi_v_phistar_ps(channel_set, 1, tau_counter) = p(1, 2);
                    
                    [r, p] = corrcoef(phi_star_partition_mis(:), phi_three_partition_phis(:));
                    phi_v_mi(channel_set, 1, tau_counter) = r(1, 2);
                    phi_v_mi_ps(channel_set, 1, tau_counter) = p(1, 2);
                    
                    [r, p] = corrcoef(phi_star_partition_mistars(:), phi_three_partition_phis(:));
                    phi_v_mistar(channel_set, 1, tau_counter) = r(1, 2);
                    phi_v_mistar_ps(channel_set, 1, tau_counter) = p(1, 2);
                    
                    % Correlate after average across trials, within conditions
                    % reshape to average within conditions
                    nTrials = size(phi_threes{nChannels_counter}.state_counters, 3);
                    % (trials*conditions x partitions) into (partitions x trials x conditions
                    phi_star_partition_phistars = reshape(phi_star_partition_phistars', [nPartitions*2, nTrials, length(conditions)]);
                    % (partitions x trials x conditions) into (partitions x conditions x trials)
                    phi_star_partition_phistars = permute(phi_star_partition_phistars, [1 3 2]);
                    % mean across trials
                    phi_star_partition_phistars = mean(phi_star_partition_phistars, 3);
                    
                    phi_star_partition_mis = reshape(phi_star_partition_mis', [nPartitions*2, nTrials, length(conditions)]);
                    phi_star_partition_mis = permute(phi_star_partition_mis, [1 3 2]);
                    phi_star_partition_mis = mean(phi_star_partition_mis, 3);
                    
                    phi_star_partition_mistars = reshape(phi_star_partition_mistars', [nPartitions*2, nTrials, length(conditions)]);
                    phi_star_partition_mistars = permute(phi_star_partition_mistars, [1 3 2]);
                    phi_star_partition_mistars = mean(phi_star_partition_mistars, 3);
                    
                    % Repeat with phi3
                    phi_three_partition_phis = reshape(phi_three_partition_phis', [nPartitions*2, nTrials, length(conditions)]);
                    phi_three_partition_phis = permute(phi_three_partition_phis, [1 3 2]);
                    phi_three_partition_phis = mean(phi_three_partition_phis, 3);
                    
                    % Correlate
                    [r, p] = corrcoef(phi_star_partition_phistars(:), phi_three_partition_phis(:));
                    trial_averaged.phi_v_phistar(channel_set, 1, tau_counter) = r(1, 2);
                    trial_averaged.phi_v_phistar_ps(channel_set, 1, tau_counter) = p(1, 2);
                    
                    [r, p] = corrcoef(phi_star_partition_mis(:), phi_three_partition_phis(:));
                    trial_averaged.phi_v_mi(channel_set, 1, tau_counter) = r(1, 2);
                    trial_averaged.phi_v_mi_ps(channel_set, 1, tau_counter) = p(1, 2);
                    
                    [r, p] = corrcoef(phi_star_partition_mistars(:), phi_three_partition_phis(:));
                    trial_averaged.phi_v_mistar(channel_set, 1, tau_counter) = r(1, 2);
                    trial_averaged.phi_v_mistar_ps(channel_set, 1, tau_counter) = p(1, 2);
                end
                
            end
            
            % Store
            correlations{nChannels_counter} = struct();
            correlations{nChannels_counter}.phi_v_phistar = phi_v_phistar;
            correlations{nChannels_counter}.phi_v_phistar_ps = phi_v_phistar_ps;
            correlations{nChannels_counter}.phi_v_mi = phi_v_mi;
            correlations{nChannels_counter}.phi_v_mi_ps = phi_v_mi_ps;
            correlations{nChannels_counter}.phi_v_mistar = phi_v_mistar;
            correlations{nChannels_counter}.phi_v_mistar_ps = phi_v_mistar_ps;
            correlations{nChannels_counter}.trial_averaged = trial_averaged;
            
        end
        
        % Store
        fly_correlations{fly_counter} = correlations;
    end
    
end

%% Concatenate flies

correlation_matrices = {...
    'phi_v_phistar',...
    'phi_v_mi',...
    'phi_v_mistar',...
    };

correlations = cell(1, length(nChannels_set));
for nChannels_counter = 1 : length(nChannels_set)
    
    % First fly
    correlations{nChannels_counter} = fly_correlations{1}{nChannels_counter};
    correlations{nChannels_counter}.nChannels = nChannels_set(nChannels_counter);
    correlations{nChannels_counter}.taus = taus;
    
    % Other flies
    for fly_counter = 2 : length(flies)
        for correlation_type = 1 : length(correlation_matrices)
            var = correlation_matrices{correlation_type};
            correlations{nChannels_counter}.(var) = cat(2, correlations{nChannels_counter}.(var), fly_correlations{fly_counter}{nChannels_counter}.(var));
            correlations{nChannels_counter}.([var '_ps']) = cat(2, correlations{nChannels_counter}.([var '_ps']), fly_correlations{fly_counter}{nChannels_counter}.([var '_ps']));
            correlations{nChannels_counter}.trial_averaged.(var) = cat(2, correlations{nChannels_counter}.trial_averaged.(var), fly_correlations{fly_counter}{nChannels_counter}.trial_averaged.(var));
            correlations{nChannels_counter}.trial_averaged.([var '_ps']) = cat(2, correlations{nChannels_counter}.trial_averaged.([var '_ps']), fly_correlations{fly_counter}{nChannels_counter}.trial_averaged.([var '_ps']));
        end
    end
    
end

%% Save

save([results_directory results_file], 'correlations');

%% Functions

function [fly_phis] = load_nChannels(fly, nChannels_set, fly_separated)
% For loading phi3 results
% Loads the relevant data for the specified nChannels for a single fly
% Concatenates nChannels together
% Assumes each fly results file has only one nChannels, i.e. length(phis)
% == 1
%
% Inputs:
%   fly - integer
%   nChannels_set - vector of nChannels which should be loaded for the fly
%   fly_separated - vector of nChannels which are separated by fly (subset
%       nChannels_set)
%
% Outputs:
%   phis - standard phi results cell array, each cell corresponds to an
%       nChannels
data_directory = 'results/preformatted_results/';

% filename format: prefix nChanneltnChannel infix fly suffix
prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels';
infix = '_phithree_allPartitionsfly';
suffix = '.mat';

fly_phis = cell(1, length(nChannels_set));

for nChannels_counter = 1 : length(nChannels_set)
    nChannels = nChannels_set(nChannels_counter);
    channel_string = [num2str(nChannels) 't' num2str(nChannels)];
    if max(ismember(fly_separated, nChannels))
        load([data_directory prefix channel_string infix num2str(fly-1) suffix]); % fly-1 is due to python's 0 indexing
    else
        load([data_directory prefix channel_string '_phithree_allPartitions' suffix]);
    end
    fly_phis{nChannels_counter} = phis{1};
end

end

%% FUNCTION

function [bipartition] = bipartition_check(partition)
% Checks if a given partition is a bipartition
%
% Inputs:
%   partition - cell array
%
% Outputs:
%   bipartition - integer; 1 if partition is bipartition, 0 otherwise

bipartition = 0;
if length(partition) == 2
    bipartition = 1;
end

end

%% Function: convert partition channel labelling from 1-indexing to 0-indexing

function [partition] = relabel(partition)
%
%
% Inputs:
%
% Outputs:

if isa(partition, 'int64')
    partition = partition_mat2cell(partition);
end

for group = 1 : length(partition)
    for channel_counter = 1 : length(partition{group})
        partition{group}(channel_counter) = partition{group}(channel_counter) + 1;
    end
end

end

function [partition_cell] = partition_mat2cell(partition)
% Converts an 2 x n partition into a cell array
%
% Inputs:
%   partition: matrix, 2 rows such that partition(1,:)--/-->(partition(2,:)
%
% Outputs:
%   partition_cell: cell array, partition{1}--/-->partition{2}

partition_cell = cell(1, size(partition, 1));

for row = 1 : size(partition, 1) % should only be 2 rows - from and to
    partition_cell{row} = partition(row, :);
end

end
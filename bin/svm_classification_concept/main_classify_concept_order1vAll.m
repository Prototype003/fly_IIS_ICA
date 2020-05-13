function [] = main_classify_concept_order1vAll(constellation_type, class_type)
% Classify using sets of concepts to match a higher order concept
% e.g. for 4th order concept ABCD, classify using 1st order A,B,C,D
%
% Classify using all 1st order concepts, and using all 2nd, 3rd, 4th order
% concepts, to check classification is completely driven by 1st orders
%
% Outputs:
%

%% Setup

% 'unpart': unpartitioned; 'part': partitioned
% 'diff': unpartitioned-partitioned; 'both': unpartitioned AND partitioned
%constellation_type = 'unpart';

%addpath('C:\Users\this_\Documents\MATLAB\Toolboxes\liblinear-2.20\windows');

addpath('../svm_classification/');

results_location = 'results/';

nChannels = 4;
flies = (1:13);
conditions = (1:2);
tau = 4;
trials = (1:8);
global_tpm = 0;
tau_bin = 0;
sample_offset = 0;

% split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f03c2tau4tauOffset0s0939t6.mat
file_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels' num2str(nChannels)];

nStates = 2^nChannels;
nConcepts = 15;

%% Load

source_dir = '../phi_3/results/';
source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels' num2str(nChannels) '_globalTPM0.mat'];
tic
disp('loading');
tmp = load([source_dir source_file]);
disp('loaded');
toc

phis = tmp.phis{1};

%% Preprocess

% Weighted average across states
big_mips = phis.big_mips;
for const_type = 1 : size(phis.big_mips, 2)
    for concept = 1 : size(phis.big_mips, 3)
        big_mips(:, const_type, concept, :, :, :, :) = permute(big_mips(:, const_type, concept, :, :, :, :), [1 4 5 6 7 2 3]) .* single(phis.state_counters);
    end
end
big_mips = permute(sum(big_mips, 1), [2 3 4 5 6 7 1]);
big_mips = big_mips ./ single(sum(phis.state_counters(:, 1, 1, 1, 1))); % State count should be constant across all parameters

% Get desired constellation type
% And specify index-vector (start of each constellation, 0-indexed as any
%   counters will be added to it, i.e. const_starts+counter)
if strcmp(constellation_type, 'diff')
    % Unpartitioned - partitioned
    big_mips = permute(...
        big_mips(1, :, :, :, :, :) - big_mips(2, :, :, :, :, :),...
        [2 3 4 5 6 1]...
        );
    const_starts = [0];
elseif strcmp(constellation_type, 'both')
    % Both unpartitioned and partitioned
    big_mips = permute(...
        cat(2, big_mips(1, :, :, :, :, :), big_mips(2, :, :, :, :, :)),...
        [2 3 4 5 6 1]...
        );
    const_starts = [0 size(phis.big_mips, 3)];
elseif strcmp(constellation_type, 'unpart')
    % Unpartitioned
    big_mips = permute(...
        big_mips(1, :, :, :, :, :),...
        [2 3 4 5 6 1]...
        );
    const_starts = [0];
elseif strcmp(constellation_type, 'part')
    % Unpartitioned
    big_mips = permute(...
        big_mips(2, :, :, :, :, :),...
        [2 3 4 5 6 1]...
        );
    const_starts = [0];
end

if strcmp(class_type, 'across')
    var_source_dim = 3; % trial dimension
else % if strcmp(class_type, 'within')
    var_source_dim = 4; % fly dimension
end

%% Classify

results_file = [num2str(nChannels) 'ch_phi3Concept_order1vAll_' constellation_type '_svm_' class_type '.mat'];

cost_powers = (-20:10:20);
costs = 2 .^ cost_powers;

% Which features to use in each SVM
order_concepts = {(1:4), (5:15)}; % hardcoded for 1st order vs 2-3-4th orders

% flies x networks x concepts+Phi x costs
cost_accuracies = zeros(size(big_mips, 2), size(big_mips, var_source_dim), length(order_concepts), length(costs));

% Create local cluster
pc = parcluster('local');

% Set JobStorageLocation to specific directory for this particular job
%pc.JobStorageLocation = strcat('matlab_pct/', getenv('SLURM_JOB_ID'));

% Start pool
%parpool(pc, 8)

% Broadcast variables
phis_p = parallel.pool.Constant(phis.phis);
big_mips_p = parallel.pool.Constant(big_mips);
const_starts_p = parallel.pool.Constant(const_starts);
order_concepts_p = parallel.pool.Constant(order_concepts);
costs_p = parallel.pool.Constant(costs);

parfor network = 1 : size(big_mips, 2)
    disp(network); tic;
    
    accuracies = zeros(size(big_mips, var_source_dim), length(order_concepts_p.Value), length(costs_p.Value));
    
    for cost_counter = 1 : length(costs_p.Value)
        cost = costs_p.Value(cost_counter);
        
        for var_instance = 1 : size(big_mips, var_source_dim)
            
            for c_order = 1 : length(order_concepts)
                
                concept_inds_tmp = order_concepts_p.Value{c_order}'; % n x 1
                concept_inds =...
                    repmat(const_starts_p.Value, [size(concept_inds_tmp, 1) 1]) +...
                    repmat(concept_inds_tmp, [1 size(const_starts_p.Value, 2)]);
                
                if strcmp(class_type, 'across')
                    
                    % Get features (concepts), average across trials
                    features = permute(...
                        big_mips_p.Value(concept_inds, network, var_instance, :, :),...
                        [4 1 5 2 3]...
                        ); % flies x features x conditions (observations x features x classes)
                    
                else % if strcmp(class_type, 'within')
                    
                    % Get features (concepts), average across trials
                    features = permute(...
                        big_mips_p.Value(concept_inds, network, :, var_instance, :),...
                        [3 1 5 2 4]...
                        ); % trials x features x conditions (observations x features x classes)
                    
                end
                
                % Classify
                results = svm_lol_liblinear_manual(features, cost);
                accuracies(var_instance, c_order, cost_counter) = results.accuracy;
                
            end
            
        end
        
    end
    
    cost_accuracies(network, :, :, :) = accuracies;
    
    toc
end


%% Save accuracies

save([results_location results_file], 'cost_accuracies', 'costs', 'cost_powers', 'nChannels', 'tau', 'order_concepts');

disp(['saved ' class_type]);

end
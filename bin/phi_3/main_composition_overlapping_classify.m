%% DESCRIPTION

%{

Plots composition Hasse graph (phi-3 is z-axis)

Average across all flies, for non-global TPM

See figures/videos from http://www.eneuro.org/content/4/5/ENEURO.0085-17.2017

%}

%% Setup

flies = (1:13);

marker_size = 100;
nChannels = 4;

addpath('../');

%% Load

tic;
%load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM0.mat');
channel_sets = double(phis{1}.channel_sets);
toc

%% Get state-weighted compositions for all parameters

composition_phis = phis{1}.big_mips;

% Weight by state occurences (multiply phi by number of times the state occurred)
for partitioned = 1 : 2
    for concept = 1 : 15
        composition_phis(:, partitioned, concept, :, :, :, :) = ...
            permute(composition_phis(:, partitioned, concept, :, :, :, :), [1 4 5 6 7 2 3]) .* ...
            double(phis{1}.state_counters);
    end
end

% Sum across states
composition_phis = permute(sum(composition_phis, 1), [2 3 4 5 6 7 1]);

% Divide by total number of states (for weighted average)
% Assumes equal number of samples for all parameters
composition_phis = composition_phis ./ sum(phis{1}.state_counters(:, 1, 1, 1, 1));

% Unpartitioned - partitioned
%composition_phis = permute(composition_phis(1, :, :, :, :, :) - composition_phis(2, :, :, :, :, :), [2 3 4 5 6 7 1]);
% Unpartitioned
composition_phis = permute(composition_phis(1, :, :, :, :, :), [2 3 4 5 6 7 1]);
% Partitioned
%composition_phis = permute(composition_phis(2, :, :, :, :, :), [2 3 4 5 6 7 1]);

%% Setup Hasse graph

% Assumes x=1 is for highest order concept, x=2 is for second highest, etc

% Find x-y space limit requirements
nConcepts = zeros(size(channel_sets, 2), 1); % Number of points needed to show all concepts
channels = max(channel_sets(:));
for concept_order = 1 : size(channel_sets, 2)
    nConcepts(concept_order) = nchoosek(channels, concept_order);
end
concepts_cum = cumsum(nConcepts);
concepts_start = concepts_cum - nConcepts + 1;

% x-y space, with padding of 1
space_min = min(channel_sets(:))-1;
space_max = max(nConcepts)+1;

% x-y space here is based on set centre and distance
x = zeros(1, sum(nConcepts));
y = zeros(size(x));
colours = zeros(size(x));
concept_counter = 1;
for concept_order = 1 : size(channel_sets, 2)
    concepts = nchoosek((1:channels), concept_order);
    for concept = 1 : length(concepts)
        x(concept_counter) = mean(concepts(concept, :)); % set center
        y(concept_counter) = channel_set_distance(concepts(concept, :)); % set distance
        concept_counter = concept_counter + 1;
    end
    colours(concepts_start(concept_order):concepts_cum(concept_order)) = size(channel_sets, 2) - concept_order + 1; % highest order at x=1
end

%% Get overlapping compositions for each fly

condition_titles = {'wake', 'anest'};

% Average across trials (result is conditions x sets x concepts x flies)
compositions = double(permute(mean(composition_phis(:, :, :, flies, :), 3), [4 5 2 1 3]));

% Find average concept phis across ALL sets which include the concept
comp_values = zeros(size(compositions, 1), length(x), size(compositions, 2));
concept_counter = 1;
concept_displacement = 0; % For skipping lower order concepts
for concept_order = 1 : size(channel_sets, 2) % 4th-order concepts aren't shared
    concepts = nchoosek((min(channel_sets(:)):max(channel_sets(:))), concept_order);
    for concept = 1 : size(concepts, 1)
        value_sum = zeros(size(compositions, 1), size(compositions, 2), 1); % values for each condition
        share_counter = 0;
        for network = 1 : size(channel_sets, 1)
            if all(ismember(concepts(concept, :), channel_sets(network, :))) % If concept is a subset of the channel set
                % Find phi of matching concept (concepts are ordered)
                network_concepts = nchoosek(channel_sets(network, :), concept_order);
                for network_concept = 1 : size(network_concepts, 1)
                    if all(ismember(concepts(concept, :), network_concepts(network_concept, :))) % If concept matches channel set concept
                        
%                         % Choose value based on max/min wake value
%                         if all(compositions(1, network, network_concept + concept_displacement) < value_sum(1))
%                             value_sum = compositions(:, network, network_concept + concept_displacement);
%                         end
                        
                        % Average values
                        value_sum = value_sum + compositions(:, :, network, network_concept + concept_displacement);
                        share_counter = share_counter + 1;
                        break;
                    end
                end
            end
        end
        
        % Selected value
%         comp_values(concept_counter, :) = value_sum;
        
        % Average values
        comp_values(:, concept_counter, :) = permute(value_sum ./ share_counter, [1 3 2]);
        
        concept_counter = concept_counter + 1;
        
    end
    concept_displacement = concept_displacement + nchoosek(size(channel_sets, 2), concept_order);
end

%% Classify across flies

addpath('C:\Users\this_\Documents\MATLAB\Toolboxes\liblinear-2.20\windows');

addpath('../svm_classification/');

cost = 1;

cost_powers = (-20:20);
costs = 2 .^ cost_powers;
cost_accuracies = zeros(size(costs));

tic;

for cost_counter = 1 : length(costs) % ~60 seconds
    cost = costs(cost_counter);
    results = svm_lol_liblinear_manual(comp_values, cost);
    cost_accuracies(cost_counter) = results.accuracy;
end

toc
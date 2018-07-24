%% Description

%{

Builds composition table (table of phis for each concept)
Currently only proof of concept, only does for one channel

%}

%% Setup

nChannels = 4;
nConcepts = 0;
for nChannels_counter = 1 : nChannels
    nConcepts = nConcepts + nchoosek(nChannels, nChannels_counter);
end
concept_list = cell(nConcepts, 1);
concept_counter = 1;
for nChannels_counter = 1 : nChannels
    mechs = nchoosek((1:nChannels), nChannels_counter);
    for mech = 1 : size(mechs, 1)
        concept_list{concept_counter} = mechs(mech, :);
        concept_counter = concept_counter + 1;
    end
end

source_dir = 'results_split/';

%% Table for wake condition

source_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM1_f01c2tau4s1036t1.mat';

load([source_dir source_file]);

nStates = length(phi.big_mips);
nMeasures = 2;
unparted_table = zeros(nStates + nMeasures, nConcepts); % n states by y mechanisms
parted_table = unparted_table;

for state = 1 : nStates
    big_mip = phi.big_mips{state};
    
    % Unpartitioned and partitioned constellations may have different number of concepts, so need separate loops
    % May want to factor out the two loops (one per constellation) into a function
    %   with inputs: constellation (unpartitioned/partitioned), state
    
    % Unpartitioned constellation
    constellation = big_mip.unpartitioned_constellation;
    concept_table_counter = 1;
    for concept_counter = 1 : length(constellation)
        concept = constellation{concept_counter};
        mech = concept.mechanism + 1;
        
        % Skip over non-existant concepts
        while length(concept_list{concept_table_counter}) ~= length(mech) || min(concept_list{concept_table_counter} == mech) ~= 1
            concept_table_counter = concept_table_counter + 1;
        end
        
        % Store phi
        unparted_table(state, concept_counter) = concept.phi;
    end
    
    % Partitioned constellation
    constellation = big_mip.partitioned_constellation;
    concept_table_counter = 1;
    for concept_counter = 1 : length(big_mip. partitioned_constellation)
        concept = constellation{concept_counter};
        mech = concept.mechanism + 1;
        
        % Skip over non-existant concepts
        while length(concept_list{concept_table_counter}) ~= length(mech) || min(concept_list{concept_table_counter} == mech) ~= 1
            concept_table_counter = concept_table_counter + 1;
        end
        
        % Store phi
        parted_table(state, concept_counter) = concept.phi;
    end
end

% Weighted sum
state = state + 1;
state_counters = permute(phi.state_counters, [2 1]);
state_counters = repmat(state_counters, [1, size(unparted_table, 2)]);
% Unpartitioned
weighted_phis = unparted_table((1:nStates), :) .* state_counters;
unparted_table(state, :) = sum(weighted_phis, 1) ./ sum(state_counters, 1);
% Partitioned
weighted_phis = parted_table((1:nStates), :) .* state_counters;
parted_table(state, :) = sum(weighted_phis, 1) ./ sum(state_counters, 1);

% Portion of existence
state = state + 1;
state_counters = permute(phi.state_counters, [2 1]);
state_counters = repmat(state_counters, [1, size(unparted_table, 2)]);
% Unpartitioned
weighted_phis = (unparted_table((1:nStates), :)>0) .* state_counters;
unparted_table(state, :) = sum(weighted_phis, 1) ./ sum(state_counters, 1);
% Partitioned
weighted_phis = (parted_table((1:nStates), :)>0) .* state_counters;
parted_table(state, :) = sum(weighted_phis, 1) ./ sum(state_counters, 1);

% (Unpartitioned - partitioned)
diff_table = unparted_table - parted_table;
% Average
weighted_phis = diff_table((1:nStates), :) .* state_counters;
dif_table(state, :) = sum(weighted_phis, 1) ./ sum(state_counters, 1);

%% Figure

concept_list_strings = concept_list;
for concept = 1 : length(concept_list)
    concept_list_strings{concept} = mat2str(concept_list{concept});
end

figure;

subplot(3, 1, 1);
imagesc(unparted_table((1:end-1), :));
title('unpartitioned');
set(gca, 'XTick', (1:15), 'XTickLabel', concept_list_strings)
set(gca, 'YTick', nStates+1, 'YTickLabel', 'mean');
ylabel('state');
c = colorbar; ylabel(c, '\phi');

subplot(3, 1, 2);
imagesc(parted_table((1:end-1), :));
title('partitioned');
set(gca, 'XTick', (1:15), 'XTickLabel', concept_list_strings)
set(gca, 'YTick', nStates+1, 'YTickLabel', 'mean');
ylabel('state');
c = colorbar; ylabel(c, '\phi');

subplot(3, 1, 3);
imagesc(diff_table((1:end-1), :));
title('partitioned - unpartitioned');
set(gca, 'XTick', (1:15), 'XTickLabel', concept_list_strings)
set(gca, 'YTick', nStates+1, 'YTickLabel', 'mean');
ylabel('state');
c = colorbar; ylabel(c, '\Delta \phi');
xlabel('concept');
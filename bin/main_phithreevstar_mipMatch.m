%% DESCRIPTION

%{
This script tests the second hypothesis: phi-3 MIP cuts are equivalent to phi-star MIPs

%}

%% SETUP

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_shareFiltered'
    ];

results_directory = 'analysis_results/';
results_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_shareFiltered'...
    '_mipComparison.mat'
    ];

%% LOAD
disp('loading');

% Phi-3
% .mips (holds the cuts) dimensions: states x sets x flies x conditions x taus
load([data_directory data_filename '_phithree.mat']);
phi_threes = phis;

% Phi-star
% .mips (holds the mips) dimensions: sets x trials x flies x conditions x taus
load([data_directory data_filename '_phistar.mat']);
phi_stars = phis;

disp('loaded')

%% Relabel phi-three MIP nodes with the channel labels
% They are labelled 0 to n-1, n is the number of channels used
% Need to map 0:n-1 to the actual channels used, given in .channel_sets

disp('relabelling');

for nChannels_counter = 1 : numel(phi_threes)
    phi_threes{nChannels_counter}.mips = relabel(phi_threes{nChannels_counter}.mips, phi_threes{nChannels_counter}.channel_sets);
end

disp('relabelled');

%% Convert phi-three cuts into partitions
% Not actually needed - due to bipartitioning only, the ends of the cuts give the partition groups

for nChannels_counter = 1 : numel(phi_threes)
    phi_threes{nChannels_counter}.mips = cuts2partitions(phi_threes{nChannels_counter}.mips);
end

%% Per trial, frequency of phi-3 MIPs matching the phi-star MIP
% Because each trial has 1 phi-star MIP and many phi-3 states (i.e. multiple MIPs), this scheme
% only allows for finding the probability of matching phi-3 MIPs for a phi-star MIP and not vice versa

matches_per_trial = cell(size(phi_threes));

% For each trial, find the ratio of matching MIPs to non matching MIPs
for nChannels_counter = 1 : numel(phi_threes)
    nChannels = phi_threes{nChannels_counter}.nChannels;
    
    % This will hold the number of matches per trial
    matches_per_trial{nChannels_counter}.matches = zeros(size(phi_stars{nChannels_counter}.mips));
    
    % This should equal the number of samples per trial (and should be constant across parameters)
    matches_per_trial{nChannels_counter}.trial_total = sum(phi_threes{nChannels_counter}.state_counters(:, 1, 1, 1, 1, 1));
    
    for tau = 1 : size(phi_stars{nChannels_counter}.mips, 5)
        for condition = 1 : size(phi_stars{nChannels_counter}.mips, 4)
            for fly = 1 : size(phi_stars{nChannels_counter}.mips, 3)
                for trial = 1 : size(phi_stars{nChannels_counter}.mips, 2)
                    for set = 1 : size(phi_stars{nChannels_counter}.mips, 1)
                        disp(['nChannels' num2str(nChannels) ' Tau' num2str(tau) ' Condition' num2str(condition) ' Fly' num2str(fly) ' trial' num2str(trial) ' set' num2str(set)]);
                        
                        % phi-star's partition
                        mip_star = phi_stars{nChannels_counter}.mips{set, trial, fly, condition, tau};
                        
                        % Determine which state(s) have the same partition
                        mip_match_indices = match_mips(mip_star, phi_threes{nChannels_counter}.mips(:, set, fly, condition, tau));
                        
                        % Count (sum; .state_counters holds the count occurrences per state) the number of matching occurrences
                        matches_per_trial{nChannels_counter}.matches(set, trial, fly, condition, tau) = sum(phi_threes{nChannels_counter}.state_counters(mip_match_indices, set, trial, fly, condition, tau));
                        
                    end
                end
            end
        end
    end
end

%% Select most common phi-three MIP as the MIP for each trial

disp('Choosing most frequent phi-3 MIPs');

phi_threes_moded = phi_threes;
for nChannels_counter = 1 : numel(phi_threes_moded)
    [phi_threes_moded{nChannels_counter}.mips, phi_threes_moded{nChannels_counter}.mip_counts] = mip_mode(phi_threes{nChannels_counter}.mips, phi_threes{nChannels_counter}.state_counters);
end

disp('Chosen');

%% Per trial, match the MIP of the most frequent state within the trial
% Because we end up having a 1:1 correspondence of phi-star to phi-3 MIPs, we can find the
% overall probability across trials of phi-star MIPs having corresponding phi-3 MIPs, and vice versa
% (probability of phi-3 MIPs having corresponding phi-star MIPs)

matches = cell(size(phi_threes));

for nChannels_counter = 1 : numel(phi_stars)
    channels_used = phi_stars{nChannels_counter}.nChannels;
    matches{nChannels_counter} = zeros(size(phi_stars{nChannels_counter}.mips)); % Equality matrix
    
    for tau = 1 : size(phi_stars{nChannels_counter}.mips, 5)
        for condition = 1 : size(phi_stars{nChannels_counter}.mips, 4)
            for fly = 1 : size(phi_stars{nChannels_counter}.mips, 3)
                for trial = 1 : size(phi_stars{nChannels_counter}.mips, 2)
                    for set = 1 : size(phi_stars{nChannels_counter}.mips, 1)
                        disp(['nChannels' num2str(channels_used) ' Tau' num2str(tau) ' Condition' num2str(condition) ' Fly' num2str(fly) ' trial' num2str(trial) ' set' num2str(set)]);
                        
                        matches{nChannels_counter}(set, trial, fly, condition, tau) = mip_equal(phi_stars{nChannels_counter}.mips{set, trial, fly, condition, tau}, phi_threes_moded{nChannels_counter}.mips{set, trial, fly, condition, tau});
                    end
                end
            end
        end
    end
    
end


%% Save the results

disp('Saving');
if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename], 'matches_per_trial', 'matches');
disp('Saved');

%% Functions: relabel phi3 channels with actual channels

function [mips] = relabel(mips, channel_sets)

% Increment channel sets to account for python's 0 indexing
channel_sets = channel_sets + 1;

for set = 1 : size(channel_sets, 1)
    
    for state = 1 : size(mips, 1)
        
        for tau = 1 : size(mips, 5)
            for condition = 1 : size(mips, 4)
                for fly = 1 : size(mips, 3)
                    partition = mips{state, set, fly, condition, tau};
                    
                    for group = 1 : length(partition)
                        partition{group} = channel_sets(set, partition{group}+1); % Increment to account for python's 0 indexing
                    end
                    
                    mips{state, set, fly, condition, tau} = partition;
                end
            end
        end
        
    end
end

end

%% Function: convert phi3 cuts into partitions

function [mips] = cuts2partitions(cuts)
% Converts phi-3 cuts to partitions (like in phi-star)
% Loop wrapper for cut2partition()
% Inputs:
%   cuts = phi_threes{n}.mips (which holds the cuts for the MIP)
% Outputs:
%   mips = cell array holding the partitions instead of the cuts (same dimensions as input cuts)

mips = cell(size(cuts));

for cut_counter = 1 : numel(cuts)
    mips{cut_counter} = cut2partition(cuts{cut_counter});
end

    function [partition] = cut2partition(cut)
        % Converts a cut into a partition
        %
        % Because cuts only give a bipartition, the partition has two groups:
        %   G1 = source nodes of the cut
        %   G2 = sink nodes of the cut
        %   e.g. for cut A > B, G1=[A] and G2=[B]
        %
        % Assumes that all nodes are referred to in the cut
        %
        % See Calculating Phi tutorial (Mayner) slide 188
        %
        % Inputs:
        %   cut = cell array with unidirectional cut of connections from elements in cut{1} to cut{2}
        % Outputs:
        %   partition = cell array with each cell holding elements which are in the same partition
        
        % Because of the bipartitioning, with all nodes referred to, the cut directly gives the partition
        
        partition = cut;
        
    end

end

%% Function: given a MIP, find states with corresponding MIPs

function [matches] = match_mips(mip, mip_set)
% Given a MIP and a set of test partitions, checks if each test partiition is equivalent to the MIP
%
% Inputs:
%   mip = cell vector holding the partitioning scheme
%   mip_set = cell vector, each cell holds the MIP of a state (see mip, above)
% Outputs:
%   matches = logical indexes of mip_set which match the mip

matches = zeros(length(mip_set), 1);
for state = 1 : length(mip_set) % For each test partition
    matches(state) = mip_equal(mip, mip_set{state});
end

matches = logical(matches);

end

%% Function: compare two MIPs

function [equal] = mip_equal(a, b)
% Compares 2 MIPs for equality
%
% Inputs:
%   a = MIP; cell array, each cell holds a vector of nodes
%   b = MIP; cell array, each cell holds a vector of nodes
%
% Outputs:
%   equal = 1 if a==b, 0 otherwise

equal = 0;

% Compare each group of MIPA to each group of MIPB
% If there is a match, we can remove the matching group from the test partition so that we don't need to look at it again
% If the B becomes empty (all groups removed), and matches is the size of the MIPA, then all groups matched
match_count = 0;
for group_a = 1 : length(a)
    for group_b = 1 : length(b)
        if length(a{group_a}) == length(b{group_b}) % Check if groups are same size
            match_group = sum(sort(a{group_a}) == sort(b{group_b})); % Check if groups have the same nodes (after sorting)
            if match_group > 0
                b{group_b} = -1; % Mark group for deletion
                match_count = match_count + 1;
            end
        end
    end
end

% Delete matched groups in the test partition (delete in reverse to avoid out-of-index errors)
for group = length(b) : -1 : 1
    if b{group} == -1
        b(group) = [];
    end
end

if match_count == length(a) && isempty(b)
    equal = 1;
end

end
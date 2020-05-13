% 'unpart': unpartitioned; 'part': partitioned
% 'diff': unpartitioned-partitioned; 'both': unpartitioned AND partitioned
constellation_type = 'unpart';

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

%%

% Get network with largest big-phi
[phi_max, max_ind] = max(phis.phis, [], 1);

% Logical matrix, 1 if network gives largest big-phi
max_status = phis.phis == repmat(phi_max, [size(phis.phis, 1) 1 1 1]);

% Get corresponding mips
big_mip_max = big_mips(repmat(permute(max_status, [5 1 2 3 4]), [size(big_mips, 1) 1 1 1 1]));
big_mip_max = reshape(big_mip_max, [size(big_mips, 1) 1 size(big_mips, 3) size(big_mips, 4) size(big_mips, 5)]);

features = permute(big_mip_max(:, 1, :, 1, :), [3 1 5 2 4]);

results = svm_lol_liblinear_manual(features, 2^-20);

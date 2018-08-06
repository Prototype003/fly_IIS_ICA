%% Description

%{

Joins 4-channel results (across all parameters)

%}

%% Setup

nChannels = 4;
global_tpm = 0;
flies = (1:13);
conditions = (1:2);
taus = [4];
tau_type = 'step'; % 'step' or 'bin'
tau_offset = 0; % Currently this script only works for a single tau_offset
trials = (1:8);

% _nChannels4_globalTPM1_f01c2tauBin4500tauOffset21s1036t1

source_dir = '../results_split/';
source_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels'...
    num2str(nChannels)...
    '_globalTPM' num2str(global_tpm)...
    ];

if strcmp(tau_type, 'step')
    tau_string = 'tau';
else % strcmp(tau_type, 'bin')
    tau_string = 'tauBin';
end

channel_sets = nchoosek((1:15), nChannels);
nStates = 2^nChannels;
dims_state_ind = [size(channel_sets, 1) length(trials), length(flies), length(conditions), length(taus)];
dims_state_dep = [nStates dims_state_ind];
nConcepts = 0; % power-set of nChannels
for concept_order = 1 : nChannels
    nConcepts = nConcepts + nchoosek(nChannels, concept_order);
end

output_file = [source_prefix(1:60) '_phithree' source_prefix(61:end) '.mat'];

%% Join split results

% NOTE: currently the actual parameters are used for indexing, so they need
%   from 1 and increment by 1
phis = cell(1);
phis{1} = struct();
phis{1}.nChannels = int8(nChannels);
phis{1}.channel_sets = int8(channel_sets);
phis{1}.taus = int8(taus);

phis{1}.phis = single(zeros(dims_state_ind));
phis{1}.big_mips = cell(dims_state_dep);
phis{1}.big_mips = single(zeros([nStates 2 nConcepts dims_state_ind]));
phis{1}.state_counters = int16(zeros(dims_state_dep));
phis{1}.state_phis = single(zeros(dims_state_dep));
phis{1}.tpms = single(zeros([nStates dims_state_dep])); % Additional dimensions because TPM is state-by-state
for fly = flies
    disp(['fly' num2str(fly)]);
    for condition = conditions
        for tau = taus
            for set_counter = 1 : size(channel_sets, 1)
                disp(['set' num2str(set_counter)]);
                for trial = trials
                    
                    % File name
                    source_file = [source_prefix...
                        '_f' sprintf('%02d', fly)...
                        'c' num2str(condition)...
                        tau_string num2str(tau)...
                        'tauOffset' num2str(tau_offset)...
                        's' sprintf('%04d', set_counter)...
                        't' num2str(trial)...
                        '.mat'...
                        ];
                    
                    try
                        % Load file
                        tmp = load([source_dir source_file]);
                    catch ME
                        disp(ME.message);
                        disp(['Failed file: ' source_file]);
                        continue; % Skip to the next file, leave the data entry structure for this entry empty
                    end
                    
                    % Place into large data structure
                    phis{1}.phis(set_counter, trial, fly, condition, tau) = single(tmp.phi.phi);
                    phis{1}.state_counters(:, set_counter, trial, fly, condition, tau) = int16(tmp.phi.state_counters);
                    phis{1}.big_mips(:, :, :, set_counter, trial, fly, condition, tau) = constellation_parse(tmp.phi.big_mips);
                    phis{1}.state_phis(:, set_counter, trial, fly, condition, tau) = single(tmp.phi.state_phis);
                    phis{1}.tpms(:, :, set_counter, trial, fly, condition, tau) = single(tmp.phi.tpm);
                    
                end
            end
        end
    end
end

%% Save

save (output_file, 'phis', '-v7.3');

%% Parse constellation

function stripped = constellation_parse(big_mips)
% Goes through big_mip and gets only the vital stuff
% Assumes that concepts are sorted!!!! This is important as details of the
% mechanism are discarded - only phi values for each concept are kept and returned
%
% Inputs:
%   big_mips: cell vector holding big_mip structs (output from pyphi.compute.big_mip)
%
% Outputs:
%   stripped: cell vector holding big_mip structs with minimal info and details
%   stripped: matrix (nMips x 2 x concepts)
%       nMips: length of 'big_mips' (for 4ch, 16)
%       2: unpartitioned (1), and partitioned (2) constellations
%       concepts: number of possible concepts (for 4ch, 15)

stripped = single(zeros(length(big_mips), 2, length(big_mips{1}.partitioned_constellation)));

for mip_counter = 1 : length(big_mips)
    
    for concept =  1 : length(big_mips{mip_counter}.unpartitioned_constellation)
        % unpartitioned_constellation
        stripped(mip_counter, 1, concept) = single(big_mips{mip_counter}.unpartitioned_constellation{concept}.phi);
        % unpartitioned_constellation
        stripped(mip_counter, 2, concept) = single(big_mips{mip_counter}.partitioned_constellation{concept}.phi);
    end
    
end

end
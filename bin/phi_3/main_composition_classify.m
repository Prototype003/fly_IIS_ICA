%% Description

%% Setup

addpath('../');

nChannels = 4;
fly = 1;
conditions = (1:2);
set = 2;
tau = 4;
trials = (1:8);
global_tpm = 0;
tau_bin = 0;
sample_offset = 0;

% split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f03c2tau4tauOffset0s0939t6.mat
file_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels' num2str(nChannels)];



%% Load compositions for each trial

for condition = conditions
    for trial = trials
        
        file_infix = [file_prefix...
            '_globalTPM' num2str(global_tpm)...
            '_f' sprintf('%02d', fly)...
            'c' num2str(condition)...
            'tau' num2str(tau)...
            'tauOffset' num2str(sample_offset)...
            's' sprintf('%04d', set)...
            't' num2str(trial)...
            '.mat'
            ];
        
        [phi, diff_table, unparted, parted] = composition_table(nChannels, [file_prefix file_infix]);
        
        
    end
end
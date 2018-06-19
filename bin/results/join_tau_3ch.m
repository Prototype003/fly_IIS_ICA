%% Description

%{

Joins 3ch phi-three 2.25s results

split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal_tau4.mat
split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal_tau8t16.mat

%}

%%

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal.mat';

%% Load
% load('split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal_tau4.mat');
% phis_tau4 = phis{1};
% 
% load('split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal_tau8t16.mat');

%% Concatenate on tau dimensions

% phis{1} already has nChannels, channel_sets, and taus [8 16]

% Update taus
phis{1}.taus = [phis_tau4.taus phis{1}.taus];

phis{1}.phi_threes = single(cat(5, phis_tau4.phi_threes, phi_threes));
phis{1}.state_counters = int16(cat(6, phis_tau4.state_counters, state_counters));
phis{1}.state_phis = single(cat(6, phis_tau4.state_phis, state_phis));
phis{1}.tpms = single(cat(7, phis_tau4.tpms, tpms));
phis{1}.mips = mips_to_int(cat(6, phis_tau4.mips, mips));

%% Save

%whos

save(results_filename, 'phis');

%% Conversions for MIPS and partitions
% Taken directly from main_parallel_phi_join.m

function [mips] = mips_to_int(mips)

for state = 1 : length(mips)
    state_mip = mips{state};
    if iscell(state_mip)
        for part = 1 : length(state_mip)
            state_mip{part} = int8(state_mip{part});
        end
    else
        state_mip = int8(state_mip);
    end
    mips{state} = state_mip;
end

end
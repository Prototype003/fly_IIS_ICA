%% Description

%{

Joins 4ch phi-three 2.25s results

split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal_tau4.mat
split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal_tau8t16.mat

%}

%%

results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal.mat';

%% Load
tau4 = load('split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal_tau4.mat');

tau_more = load('split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal_tau8t16.mat');

%% Concatenate on tau dimensions

phis = cell(1);
phis{1}.nChannels = 4;
phis{1}.channel_sets = tau4.phis{1}.channel_sets;
phis{1}.taus = [tau4.taus tau_more.taus];

flies = (1:size(tau_more.phi_threes, 3));
phis{1}.phi_threes = single(cat(5, tau4.phi_threes(:, :, flies, :), tau_more.phi_threes));
phis{1}.state_counters = int16(cat(6, tau4.state_counters(:, :, :, flies, :), tau_more.state_counters));
phis{1}.state_phis = single(cat(6, tau4.state_phis(:, :, :, flies, :), tau_more.state_phis));
phis{1}.tpms = single(cat(7, tau4.tpms(:, :, :, :, flies, :), tau_more.tpms));
%phis{1}.mips = mips_to_int(cat(6, tau4.mips(:, :, :, flies, :), tau_more.mips));

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
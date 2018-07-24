%% Description

%{

Computes phi values

%}

%% Settings

set_sizes = (2:4);
phi_type = 'star';

flies = (1:13);
conditions = (1:2);
trials = (1:8);
taus = [4 8 16];

channels = (1:15);

results_dir = 'results_split/';

%% Compute and save

phis = cell(length(set_sizes), 1);
for nChannels_counter = 1 : length(set_sizes)
    nChannels = set_sizes(nChannels_counter);
    channel_sets = nchoosek(channels, nChannels);
    sets = (1:size(channel_sets, 1));
    
    results_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels' num2str(nChannels) '_phi' phi_type '_nonGlobal.mat'];
    phis{nChannels_counter}.nChannels = nChannels;
    phis{nChannels_counter}.channel_sets = channel_sets;
    phis{nChannels_counter}.taus = taus;
    phis{nChannels_counter}.phis = zeros(length(sets), length(trials), length(flies), length(conditions), length(taus));
    phis{nChannels_counter}.mips = cell(length(sets), length(trials), length(flies), length(conditions), length(taus));
    for fly = 1 : length(flies)
        for condition = 1 : length(conditions)
            for tau_counter = 1 : length(taus)
                tau = taus(tau_counter);
                for set_counter = 1 : length(sets)
                    for trial = 1 : length(trials)
                        disp(['nChannels:' num2str(nChannels) ' fly:' num2str(fly) ' cond:' num2str(condition) ' tau:' num2str(tau) ' set:' num2str(set_counter) ' trial:' num2str(trial)]);
                        [mip, phi] = compute_discrete_nonGlobal(phi_type, fly, condition, trial, nChannels, set_counter, tau, 0);
                        phis{1}.phis(set_counter, trial, fly, condition, tau_counter) = phi;
                        phis{1}.mips{set_counter, trial, fly, condition, tau_counter} = mip;
                    end
                end
            end
        end
    end
    
    save([results_dir results_file], 'phis');
end
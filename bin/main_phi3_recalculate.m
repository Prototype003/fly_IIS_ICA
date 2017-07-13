%% DESCRIPTION

%{
This script recalculates phis at each trials from the state-phis and state-counters
Sepearate fly results should have already been joined together
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
    '_phithree.mat'
    ];

results_directory = 'results/';
results_filename = data_filename; % As the existing data source has wrong phi-values, it will be overwritten

%% Load base results

disp('loading');
load([data_directory data_filename]);
disp('loaded');

%% Recalculate phis per trial

for nChannels_counter = 1 : length(phis)
    for tau = 1 : size(phis{nChannels_counter}.phi_threes, 5)
        for condition = 1 : size(phis{nChannels_counter}.phi_threes, 4)
            for fly = 1 : size(phis{nChannels_counter}.phi_threes, 3)
                for trial = 1 : size(phis{nChannels_counter}.phi_threes, 2)
                    disp(['tau' num2str(tau) ' condition' num2str(condition) ' fly' num2str(fly) ' trial' num2str(trial)]);
                    for set = 1 : size(phis{nChannels_counter}.phi_threes, 1)
                        phis{nChannels_counter}.phi_threes(set, trial, fly, condition, tau) =...
                            sum(phis{nChannels_counter}.state_phis(:, set, fly, condition, tau) .* phis{nChannels_counter}.state_counters(:, set, trial, fly, condition, tau)) /...
                            sum(phis{nChannels_counter}.state_counters(:, set, trial, fly, condition, tau)); % This should be constant (=number of samples = 2250)
                    end
                end
            end
        end
    end
end

%% Save (overwrite old file with wrong phi values)
% Save directory is the load directory, so it definitely exists
disp('Saving');
save([results_directory results_filename], 'phis');
disp('Saved');
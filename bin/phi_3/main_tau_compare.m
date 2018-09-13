%% Description

%{
Compares phi values for tau values

for N flies, Y channel sets, etc

%}

%% Setup

nChannels = 4;
flies = (1:13);
conditions = (1:2);
taus = 2.^(0:9);
%taus = 2.^(0:6);

nchoosek(15, nChannels);
sets = floor(linspace(1, nchoosek(15, nChannels), 10));

trials = (1);
binOffsets = 0;

global_tpm = 1;
tau_string = 'tau'; % 'tauBin' or 'tau'

% 'tmp/' = results without any correction
% 'tmp2/' = results with TPM 'normalisation' by penalising transition
%   probabilities with less than n observations
% 'tmp3/' = results with 'fairer' TPMs, built by offsetting samples before
%   binning into coarse grained samples
source_dir = 'tmp/';
source_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_nChannels' num2str(nChannels)...
    '_globalTPM' num2str(global_tpm)...
    ];

%% Load phis

% % sets x flies x conditions x taus
% phi_values = zeros(length(sets), length(flies), length(conditions), length(taus));
% 
% for set_counter = 1 : length(sets)
%     network = sets(set_counter);
%     disp(network);
%     for fly = flies
%         for condition = conditions
%             for tau_counter = 1 : length(taus)
%                 for trial = trials
%                     
%                     % File name (for tmp3/)
%                     source_file = [source_prefix...
%                         '_f' sprintf('%02d', fly)...
%                         'c' num2str(condition)...
%                         tau_string num2str(taus(tau_counter))...
%                         'binOffsets' num2str(binOffsets)...
%                         's' sprintf('%04d', network)...
%                         't' num2str(trial)...
%                         '.mat'...
%                         ];
%                     
% %                     % File name (for tmp/ and tmp2)
% %                     source_file = [source_prefix...
% %                         '_f' sprintf('%02d', fly)...
% %                         'c' num2str(condition)...
% %                         tau_string num2str(taus(tau_counter))...
% %                         'tauOffset' num2str(0)...
% %                         's' sprintf('%04d', network)...
% %                         't' num2str(trial)...
% %                         '.mat'...
% %                         ];
%                     
%                     try
%                         % Load file
%                         tmp = load([source_dir source_file]);
%                     catch ME
%                         disp(ME.message);
%                         disp(['Failed file: ' source_file]);
%                         continue; % Skip to the next file, leave the data entry structure for this entry empty
%                     end
%                     
%                     phi_values(set_counter, fly, condition, tau_counter) = tmp.phi.phi;
%                     
%                 end
%             end
%         end
%     end
% end
% 
% %% Save (so we don't need to keep loading)
% 
% save([source_dir tau_string '.mat'], 'phi_values');

%% Load the saved one instead of all the separate files

load([source_dir tau_string '.mat']);

%% Plot phi averaged across sets, flies

%phi_values = phi_values(:, :, :, 1:length(taus));

flies = (1:13);

% figure;
% imagesc(squeeze(mean(phi_values(:, flies, 1, :), 2)));
% %imagesc((squeeze(mean(phi_values(:, flies, 1, :), 2)) - squeeze(mean(phi_values(:, flies, 2, :), 2))) ./ squeeze(mean(phi_values(:, flies, 1, :), 2)));
% 

tau_log = log2(taus);

figure;

% Plot for each individual channel set
for set_counter = 1 : length(sets)
    
    subplot(2, 3, 1);
    plot(tau_log, squeeze(mean(phi_values(set_counter, flies, 1, :), 2))); hold on;
    title('wake');
    ylabel('\Phi');
    xlabel('log2(\tau)');
    xlim([tau_log(1) tau_log(end)]);
    ax = gca;
    ax.YAxis.Exponent = 0;
    yl = ylim;
    
    subplot(2, 3, 2);
    plot(tau_log, squeeze(mean(phi_values(set_counter, flies, 2, :), 2))); hold on;
    title('anest');
    xlim([tau_log(1) tau_log(end)]);
    ylim(yl);
    ax = gca;
    ax.YAxis.Exponent = 0;
    
    subplot(2, 3, 3);
    plot(tau_log, squeeze(mean(phi_values(set_counter, flies, 1, :), 2)) - squeeze(mean(phi_values(set_counter, flies, 2, :), 2))); hold on;
    title('wake-anest');
    xlim([tau_log(1) tau_log(end)]);
    ax = gca;
    ax.YAxis.Exponent = 0;
end

% Mean and std across the channel sets
subplot(2, 3, 4);
errorbar(tau_log, squeeze(mean(mean(phi_values(:, flies, 1, :), 2), 1)), squeeze(std(mean(phi_values(:, flies, 1, :), 2), [], 1))); hold on;
title('wake');
ylabel('\Phi');
xlabel('log2(\tau)');
xlim([tau_log(1)-1 tau_log(end)+1]);
ax = gca;
ax.YAxis.Exponent = 0;
yl = ylim;

subplot(2, 3, 5);
errorbar(tau_log, squeeze(mean(mean(phi_values(:, flies, 2, :), 2), 1)), squeeze(std(mean(phi_values(:, flies, 2, :), 2), [], 1))); hold on;
title('anest');
xlim([tau_log(1)-1 tau_log(end)+1]);
ax = gca;
ax.YAxis.Exponent = 0;
ylim(yl);

subplot(2, 3, 6);
errorbar(tau_log, squeeze(mean(mean(phi_values(:, flies, 1, :)-phi_values(:, flies, 2, :), 2), 1)), squeeze(std(mean(phi_values(:, flies, 1, :)-phi_values(:, flies, 2, :), 2), [], 1))); hold on;
title('wake-anest');
xlim([tau_log(1)-1 tau_log(end)+1]);
ax = gca;
ax.YAxis.Exponent = 0;

%legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10');
%% With variance across 2250ms trials

phi_values = zeros(length(taus), length(trials));
for tau_counter = 1 : length(taus)
    tau = taus(tau_counter);
    
    for trial_counter = 1 : length(trials)
        trial = trials(trial_counter);
        
        source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f01c1tauBin' num2str(tau) 's1036t' num2str(trial) '.mat'];
        
        load([source_dir source_file]);
        
        phi_values(tau_counter, trial) = phi.phi;
    end
end

phi_mean = mean(phi_values, 2);
phi_std = std(phi_values, [], 2);

figure;
errorbar(taus, phi_mean, phi_std);

figure;
plot(taus, (phi_std ./ phi_mean) ./ taus');
hold on;
plot(taus, phi_std ./ taus');

%% With variance across offsets

% phi_values = zeros(length(taus), nOffsets);
% for tau_counter = 1 : length(taus)
%     tau = taus(tau_counter);
%     
%     for offset = 1 : nOffsets
%         source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM1_f01c1tauBin' num2str(tau) 'tauOffset' num2str(offset-1) 's1036t1.mat'];
%         
%         load([source_dir source_file]);
%         
%         phi_values(tau_counter, offset) = phi.phi;
%     end
% end
% 
% phi_mean = mean(phi_values, 2);
% phi_std = std(phi_values, [], 2);
% 
% figure;
% plot(taus, (phi_std ./ phi_mean) ./ taus');
% hold on;
% plot(taus, phi_std ./ tau_bins');

%% Plots mean across variance source (across trials or offsets)

figure;
errorbar(taus, phi_mean, phi_std);
xlabel('tau bin size (ms)');
ylabel('\phi');
title(['mean of ' num2str(nOffsets) ' sample offsets']);

%% Plots every variance source (every trial or every offset)

figure;
for offset = 1 : length(trials) %nOffsets
    plot(taus, phi_values(:, offset)); hold on;
end

%% Classification performance (across trials)

addpath('../');

accuracies = zeros(length(taus), 1);

% Get data matrix
values = zeros(length(trials), 2, length(taus));
for tau_counter = 1 : length(taus)
    tau = taus(tau_counter);
    
    for trial_counter = 1 : length(trials)
        trial = trials(trial_counter);
        
        % Load awake
        source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f01c1tauBin' num2str(tau) 's1036t' num2str(trial) '.mat'];
        load([source_dir source_file]);
        values(trial_counter, 1, tau_counter) = phi.phi;
        
        % Load anest
        source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f01c2tauBin' num2str(tau) 's1036t' num2str(trial) '.mat'];
        load([source_dir source_file]);
        values(trial_counter, 2, tau_counter) = phi.phi;
    end
    
    % Classify
    [accuracies(tau_counter), leave_out, correctness] = classify_lol(values(:, :, tau_counter));
end

figure;
plot((1:length(taus)), accuracies);
set(gca, 'XTick', (1:length(taus)), 'XTickLabel', taus);
xlim([1, length(taus)]);
xlabel('tau bin size');
ylabel('classification accuracy');
title('1 fly, 1 channel set, 2.25s TPM; classification across 2x8 trials');

%% Classification performance (across offsets)

addpath('../');

accuracies = zeros(length(taus), 1);

% Get data matrix
values = zeros(nOffsets, 2, length(taus));
for tau_counter = 1 : length(taus)
    tau = taus(tau_counter);
    
    for offset = 1 : nOffsets
        % Load awake
        source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM1_f01c1tauBin' num2str(tau) 'tauOffset' num2str(offset-1) 's1036t1.mat'];
        load([source_dir source_file]);
        values(offset, 1, tau_counter) = phi.phi;
        
        % Load anest
        source_file = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM1_f01c2tauBin' num2str(tau) 'tauOffset' num2str(offset-1) 's1036t1.mat'];
        load([source_dir source_file]);
        values(offset, 2, tau_counter) = phi.phi;
    end
    
    % Classify
    [accuracies(tau_counter), leave_out, correctness] = classify_lol(values(:, :, tau_counter));
end

figure;
plot((1:length(taus)), accuracies);
set(gca, 'XTick', (1:length(taus)), 'XTickLabel', taus);
xlim([1, length(taus)]);
xlabel('tau bin size');
ylabel('classification accuracy');
title('1 fly, 1 channel set, 18s TPM; classification across 25 offsets');
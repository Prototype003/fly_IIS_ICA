%% Description

%{
Compares phi values for binned tau values

for fly 1:
    2ch : 61
    4ch : 1036
%}

%% Setup

taus = [1 2 3 4 8 12 16 24 32 48 64 128 256 512];
%taus = (1:32);

nOffsets = 25; % 4 or 25

trials = (1:8);

source_dir = 'results_split/';

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
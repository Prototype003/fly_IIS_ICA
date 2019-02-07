%% Load

results_directory = 'workspace_results/';
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_power_classification_across1.mat';

load([results_directory results_filename]);

nChannels = size(accuracies, 2);
nFlies = size(accuracies, 3);

%% Plot average

figure;
errorbar(...
    frequencies,...
    mean(mean(accuracies, 2), 3),... % Average across channel pairs, then the flies
    0.5*std(mean(accuracies, 2), [], 3) / sqrt(nFlies)); % standard error across flies
title('Classification (anest/awake) using power');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass 45 100]);

% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);

%% Plot for 1 fly

fly = 3;
figure;
for channel = 1 : nChannels
    plot(frequencies, accuracies(:, channel, fly)); hold on;
end
title('Classification (anest/awake) using power');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass 0 100]);

% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);
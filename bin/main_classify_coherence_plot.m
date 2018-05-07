%% Load

results_directory = 'workspace_results/';
results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_coherence_classification.mat';

load([results_directory results_filename]);

nNetworks = size(accuracies, 2);
nFlies = size(accuracies, 3);

%% Plot average

figure;
errorbar(...
    frequencies,...
    mean(mean(accuracies, 2), 3),... % Average across channel pairs, then the flies
    0.5*std(mean(accuracies, 2), [], 3) / sqrt(nFlies)); % Standard error across flies
title('Classification (anest/awake) using power');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass 45 100]);

% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);

%% Plot for 1 fly

fly = 11;
figure;
for network = 1 : nNetworks
    plot(frequencies, accuracies(:, network, fly)); hold on;
end
title('Classification (anest/awake) using power');
xlabel('frequency (Hz)');
ylabel('% correct');
axis([chronux_params.fpass 0 100]);

% Chance level
line(chronux_params.fpass, [100/nConditions 100/nConditions]);
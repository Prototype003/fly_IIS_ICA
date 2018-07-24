%% Description

%{

Plots channel positions of selected channel channels used in 5ch phi-three

%}

%% Setup

nChannels = 5;

bin_location = '../';

source_directory = ['results_phiChain/'];
source_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels5_globalTPM0_chainedPhi.mat';

%% Load

load([source_directory source_file]);

%% Plot channel sets per fly

% Channel sets should be identical across trials

channel_selection = zeros(size(phi_joined.channel_sets, 3), 15);

for fly = 1 : size(phi_joined.channel_sets, 3)
    channel_selection(fly, phi_joined.channel_sets(:, 1, fly)) = 1;
end

figure;
subplot(1, 4, [1 2]);
imagesc(channel_selection);

% Plot phi values
subplot(1, 4, 3);
values = permute(mean(phi_joined.phis, 1), [2 3 1]); % flies x conditions
plot(values(:, 1), (1:size(values, 1)), 'r'); hold on;

subplot(1, 4, 3);
values = permute(mean(phi_joined.phis, 1), [2 3 1]); % flies x conditions
plot(values(:, 2), (1:size(values, 1)), 'b');

ylabel('fly'); xlabel('phi3'); legend('w', 'a');

subplot(1, 4, 4);
values = phi_joined.accuracies;
plot(values, (1:size(phi_joined.phis, 2)));

ylabel('fly'); xlabel('class. %');
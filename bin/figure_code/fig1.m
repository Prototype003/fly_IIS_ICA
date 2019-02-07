%% Settings

fly = 1;
channels = [5 6]; %[5 6];
trial = 1;
samples = (101:120); %(101:110);
condition = 1;
tau = 1;

addpath('../');

%% Load
data_directory = '../workspace_results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

disp('Loading fly data');
loaded_data = load([data_directory data_file '.mat']);
fly_data = loaded_data.fly_data; % Reminder: dimensions are: (samples x channels x trials x flies x conditions)
disp('Fly data loaded')

%% Preprocess

% Binarise fly_data
%fly_data_binarised = binarise_global_median(fly_data);
n_values = 2;

% Get relevant data
raw_data = fly_data(samples, channels, trial, fly, condition);

% Each channel is binarised based on its median value
middle = median(raw_data, 1);
middle_mat = repmat(middle, [size(raw_data, 1), 1]);
binarised_data = raw_data > middle_mat;
channel_data = binarised_data; %fly_data_binarised(samples, channels, trial, fly, condition);

% Actual TPM
tpm = build_tpm(channel_data, 1, n_values);

% Independent TPM
[tpm_ind, ind_a, ind_b] = build_tpm_independent(channel_data, 1, n_values);

%% Create figure

figure;
set(gcf, 'Position', [0 0 2100/2 500]);
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
set(gcf, 'RendererMode', 'manual');
set(gcf, 'Renderer', 'painters');

text_labels = 'afcd';

%% Plot data

% Color plot
% data_plot = subplot(5, 6, (1:3));
% imagesc(raw_data'); cbar = colorbar; axis_defaults(gca);
% set(gca, 'YTick', [1 2 3], 'XTickLabel', '');
% xlabel(cbar, 'V');%'\muV');
% ylabel('channel');
% colormap(data_plot, 'jet');

% Line plot
line_colour = [0.85 0.325 0.098];
nLines = size(raw_data, 2);
line_offsets = (1:nLines)*(max(raw_data(:)) - min(raw_data(:)));
medians = zeros(size(line_offsets));
data_plot = subplot(5, 6, (1:3));
for line_counter = 1 : nLines
    line([1 size(raw_data, 1)], repmat(median(raw_data(:, line_counter))-line_offsets(line_counter), [1 2]), 'Color', [0 0 0 0.5], 'LineWidth', 1);
    medians(line_counter) = median(raw_data(:, line_counter)) - line_offsets(line_counter);
    hold on;
    plot(raw_data(:, line_counter)-line_offsets(line_counter), '-o', 'Color', line_colour, 'MarkerSize', 2, 'MarkerFaceColor', line_colour, 'LineWidth', 1); hold on;
end
line([20 20], [mean(medians)-12.5 mean(medians)+12.5], 'Color', 'k'); % Draw scale line
text(20.5, mean(medians), '25\muV');
axis_defaults(gca);
set(gca, 'YTick', sort(medians), 'YTickLabel', (nLines:-1:1), 'XTick', [], 'XColor', [1 1 1 0]);
xlim([0.5 size(raw_data, 1)+0.5]);
ylabel('channel');

fontsize = 11;
textbox_width = 0.03;
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(1), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');


binarised_plot = subplot(5, 6, (7:9));
imagesc(channel_data'); axis_defaults(gca);
set(gca, 'YTick', [1 2]);
cbar = colorbar('location', 'eastoutside'); caxis([0 1]);
set(cbar, 'YTick', [0.25 0.75], 'YTickLabel', {'off (0)', 'on (1)'});
plot_pos = get(gca, 'Position'); pos_ref = get(data_plot, 'Position');
plot_pos(3:4) = pos_ref(3:4); set(gca, 'Position', plot_pos);
xlabel('time sample'); ylabel('channel');
colormap(binarised_plot, [0 0 0; 1 1 1]);

%% Plot TPMs

states_ind = {'0', '1'};
states = {'00', '01', '10', '11'};

plim = [0 1]; % Probability colourbar

% Empirical TPM
tpm_plot = subplot(5, 6, [4.5 5.2 10.5 11.2]);
imagesc(tpm, plim);
plot_pos = get(gca, 'Position'); axis square;
cbar = colorbar('location', 'eastoutside'); xlabel(cbar, 'P'); set(cbar, 'XTick', [0 0.5 1]); axis_defaults(gca);
set(gca, 'Position', plot_pos);
cbar_pos = get(cbar, 'Position'); pos_ref = get(gca, 'Position');
cbar_pos(4) = pos_ref(4); set(cbar, 'Position', cbar_pos);
set(gca, 'XTick', [1 2 3 4], 'YTick', [1 2 3 4], 'YTickLabel', states, 'XTickLabel', states);
xlabel('state t+1'); ylabel('state t');
%colormap(tpm_plot, 'gray');
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(2), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');

% Independent channel TPMs
ind_a_plot = subplot(5, 6, [13 13.1]+6);
imagesc(ind_a, plim); cbar = colorbar; cbar.Visible = 'off'; axis_defaults(gca); axis square;
set(gca, 'XTick', [1 2], 'YTick', [1 2], 'YTickLabel', states_ind, 'XTickLabel', states_ind);
xlabel('t+1'); ylabel('t');
title('ch1');
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1)-0.01 ax_pos(2)+ax_pos(4)+0.01 textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(3), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');

%colormap(ind_a_plot, 'gray');
ind_b_plot = subplot(5, 6, [15 15.1]+6);
imagesc(ind_b, plim); cbar = colorbar; cbar.Visible = 'off'; axis_defaults(gca); axis square;
set(gca, 'XTick', [1 2], 'YTick', [1 2], 'YTickLabel', states_ind, 'XTickLabel', states_ind);
xlabel('t+1'); ylabel('t');
title('ch2');
%colormap(ind_b_plot, 'gray');

% Independent TPM (product of indepent channel TPMs)
tpm_ind_plot = subplot(5, 6, [4.5 5.2 10.5 11.2]+18);
imagesc(tpm_ind, plim);
tmp_pos = get(gca, 'Position'); axis square;
tmp_pos(3:4) = plot_pos(3:4); set(gca, 'Position', tmp_pos);
plot_pos = get(gca, 'Position');
cbar = colorbar('location', 'eastoutside'); xlabel(cbar, 'P'); set(cbar, 'XTick', [0 0.5 1]); axis_defaults(gca); set(gca, 'Position', plot_pos);
cbar_pos = get(cbar, 'Position'); pos_ref = get(gca, 'Position');
cbar_pos(4) = pos_ref(4); set(cbar, 'Position', cbar_pos);
set(gca, 'XTick', [1 2 3 4], 'YTick', [1 2 3 4], 'YTickLabel', states, 'XTickLabel', states);
xlabel('state t+1'); ylabel('state t');
%colormap(tpm_ind_plot, 'gray');
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(4), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');


%% Print figure

% figure_name = 'fig1b';
% 
% print(figure_name, '-dsvg', '-painters'); % SVG
% print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
% print(figure_name, '-dpng'); % PNG

%% Show probability distributions for a given state

state = 4;

figure;
set(gcf, 'color', 'w');

subplot(1, 2, 1);
bar(tpm(state, :), 'k');
bar([1/4 1/4 1/4 1/4], 'k');
bar([0 1 0 0], 'k');
ylim([0 1]);
xlim([0.5 4.5]);
set(gca, 'XTick', (1:4), 'XTickLabel', {'00','01','10','11'});
set(gca, 'YTick', [0 0.5 1]);
xlabel('state_{t+1}');
ylabel('P(state)', 'Rotation', 90);
box off;

subplot(1, 2, 2);
bar(tpm_ind(state, :), 'k');
ylim([0 1]);
xlim([0.5 4.5]);
set(gca, 'XTick', (1:4), 'XTickLabel', {'00','01','10','11'});
set(gca, 'YTick', [0 0.5 1]);
xlabel('state_{t+1}');
ylabel('P(state)', 'Rotation', 90);
box off;

%% Function: build independent TPM from 2 channel data
function [tpm, a, b] = build_tpm_independent(channel_data, tau, n_values)
% For the 2 channel scenario
% Multiplies two (independent) TPMs together using Kronecker Tensor
% multiplication (kron())
% Note that the output is NOT in the required LOLI format for PyPhi

% Build TPM for single channels
a = build_tpm(channel_data(:, 1), tau, n_values);
b = build_tpm(channel_data(:, 2), tau, n_values);

% Multiply single-channel TPMs
tpm = kron(a, b);

end
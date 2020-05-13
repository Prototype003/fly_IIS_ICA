%% Description

%{

Create methods figure

a) fly recording scheme to raw LFPs and binarised time course
    (same as fig1.m in ../figure_code/)
b) full SBN TPM
c) independent SBN TPM
d) establish axis (x=concept order; y=phi)
    rotated view for unpartitioned, partitioned, unpart-part composition
%}

%% Settings

fly = 1;
channels = [2 4]; %[1 2]; %[5 6];
trial = 1;
samples = (1:20); %(101:120); %(1:20);
condition = 1;
tau = 1;

addpath('../');

%% Load raw data

data_directory = '../workspace_results/';
data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim';

disp('Loading fly data');
loaded_data = load([data_directory data_file '.mat']);
fly_data = loaded_data.fly_data; % Reminder: dimensions are: (samples x channels x trials x flies x conditions)
disp('Fly data loaded')

%% Load TPM and mechanism-purview values

sample_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels2_globalTPM0_f01c1tau4s0016t1_example.mat';

load(sample_file);

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

%% Build TPM for sample range

% Full
tpm_sbs = build_tpm(channel_data, 1, 2);
tpm_sbn = tpm_sbs2sbn(tpm_sbs);

% Split
sbs_0 = build_tpm(channel_data(:, 1), 1, 2);
sbs_1 = build_tpm(channel_data(:, 2), 1, 2);
tpm_sbn_0 = tpm_sbs2sbn(sbs_0);
tpm_sbn_1 = tpm_sbs2sbn(sbs_1);
tpm_sbn_ind = tpm_sbn_product(cat(2, tpm_sbn_0, tpm_sbn_1));

% figure;
% subplot(1, 2, 1); imagesc(tpm_sbn); colorbar;
% subplot(1, 2, 2); imagesc(tpm_sbn_ind); colorbar;

%% Create figure

fig = figure;

text_labels = 'afcd';

set(gcf, 'Position', [0 0 2100/2 500]);
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
set(gcf, 'RendererMode', 'manual');
set(gcf, 'Renderer', 'painters');

%% Plot data

% Line plot
line_colour = [0.85 0.325 0.098];
nLines = size(raw_data, 2);
line_offsets = [0 0.5]*(max(raw_data(:)) - min(raw_data(:)));
medians = zeros(size(line_offsets));
data_plot = subplot(5, 8, (1:4));
for line_counter = 1 : nLines
    line([1 size(raw_data, 1)], repmat(median(raw_data(:, line_counter))-line_offsets(line_counter), [1 2]), 'Color', [0 0 0 0.5], 'LineWidth', 1);
    medians(line_counter) = median(raw_data(:, line_counter)) - line_offsets(line_counter);
    hold on;
    plot(raw_data(:, line_counter)-line_offsets(line_counter), '-o', 'Color', line_colour, 'MarkerSize', 2, 'MarkerFaceColor', line_colour, 'LineWidth', 1); hold on;
end
line([20 20], [mean(medians)-12.5 mean(medians)+12.5], 'Color', 'k'); % Draw scale line
text(20.5, mean(medians), '25\muV');
axis_defaults(gca);
set(gca, 'YTick', sort(medians), 'YTickLabel', fliplr({'A', 'B'}), 'XTick', [], 'XColor', [1 1 1 0]);
xlim([0.5 size(raw_data, 1)+0.5]);
ylabel('channel');

fontsize = 11;
textbox_width = 0.03;
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(1), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');


binarised_plot = subplot(5, 8, (1:4)+8);
imagesc(channel_data'); axis_defaults(gca);
set(gca, 'YTick', [1 2]);
cbar = colorbar('location', 'eastoutside'); caxis([0 1]);
set(cbar, 'YTick', [0.25 0.75], 'YTickLabel', {'off (0)', 'on (1)'});
plot_pos = get(gca, 'Position'); pos_ref = get(data_plot, 'Position');
plot_pos(3:4) = pos_ref(3:4); set(gca, 'Position', plot_pos);
xlabel('time sample'); ylabel('channel');
set(gca, 'YTickLabel', {'A', 'B'});
colormap(binarised_plot, [0 0 0; 1 1 1]);

%% Plot TPMs

states_ind = {'0', '1'};
states = {'00', '10', '01', '11'};

dpoints = 2; % decimal points

% Concatenate TPMs to find min and max
inds_joined = [tpm_sbn_0 tpm_sbn_1];
inds_joined = [inds_joined; tpm_sbn_ind; tpm_sbn];
plim = [min(inds_joined(:)) max(inds_joined(:))]; % Probability colourbar

% Empirical TPM
tpm_plot = subplot(5, 8, [6 14]);
imagesc(tpm_sbn, plim);
plot_pos = get(gca, 'Position'); %axis square;
cbar = colorbar('location', 'eastoutside'); xlabel(cbar, 'P'); set(cbar, 'XTick', [min(plim) 0.5 max(plim)]); set(cbar, 'XTickLabel', round([min(plim) 0.5 max(plim)], dpoints)); axis_defaults(gca);
set(gca, 'Position', plot_pos);
cbar_pos = get(cbar, 'Position'); pos_ref = get(gca, 'Position');
cbar_pos(4) = pos_ref(4); set(cbar, 'Position', cbar_pos);
set(gca, 'xaxisLocation', 'top');
set(gca, 'XTick', [1 2 3 4], 'YTick', [1 2 3 4], 'YTickLabel', states, 'XTickLabel', {'A=1','B=1'});
ylabel('state t');
%colormap(tpm_plot, 'gray');
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(2), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');

% Independent channel TPMs
ind_a_plot = subplot(5, 8, [2 2.05]+24);
imagesc(tpm_sbn_0, plim); cbar = colorbar; cbar.Visible = 'off'; axis_defaults(gca); %axis square;
set(gca, 'xaxisLocation', 'top');
set(gca, 'XTick', [1 2], 'YTick', [1 2], 'YTickLabel', states_ind, 'XTickLabel', {'A=1'});
ylabel('t');
title('ch1');
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1)-0.01 ax_pos(2)+ax_pos(4)+0.01 textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(3), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');

%colormap(ind_a_plot, 'gray');
ind_b_plot = subplot(5, 8, [4 4.05]+24);
imagesc(tpm_sbn_1, plim); cbar = colorbar; cbar.Visible = 'off'; axis_defaults(gca); %axis square;
set(gca, 'xaxisLocation', 'top');
set(gca, 'XTick', [1 2], 'YTick', [1 2], 'YTickLabel', states_ind, 'XTickLabel', {'B=1'});
ylabel('t');
title('ch2');
%colormap(ind_b_plot, 'gray');

% Independent TPM (product of indepent channel TPMs)
tpm_plot = subplot(5, 8, [6 14]+24);
imagesc(tpm_sbn_ind, plim);
tmp_pos = get(gca, 'Position'); %axis square;
tmp_pos(3:4) = plot_pos(3:4); set(gca, 'Position', tmp_pos);
plot_pos = get(gca, 'Position');
cbar = colorbar('location', 'eastoutside'); xlabel(cbar, 'P'); set(cbar, 'XTick', [min(plim) 0.5 max(plim)]); set(cbar, 'XTickLabel', round([min(plim) 0.5 max(plim)], dpoints)); axis_defaults(gca); set(gca, 'Position', plot_pos);
cbar_pos = get(cbar, 'Position'); pos_ref = get(gca, 'Position');
cbar_pos(4) = pos_ref(4); set(cbar, 'Position', cbar_pos);
set(gca, 'xaxisLocation', 'top');
set(gca, 'XTick', [1 2 3 4], 'YTick', [1 2 3 4], 'YTickLabel', states, 'XTickLabel', {'A=1','B=1'});
ylabel('state t');
%colormap(tpm_ind_plot, 'gray');
ax_pos = get(gca, 'Position');
axes('Visible', 'off', 'Position', [ax_pos(1) ax_pos(2)+ax_pos(4) textbox_width textbox_width]); axis_defaults(gca);
text(0, 0.2, text_labels(4), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize, 'FontWeight', 'bold');

%% Print

figure_name = 'figures/fig1a';

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG

%% Plot TPMs (as bar graphs)

states_ind = {'0', '1'};
states = {'00', '10', '01', '11'};

dpoints = 2; % decimal points

% Concatenate TPMs to find min and max
inds_joined = [tpm_sbn_0 tpm_sbn_1];
inds_joined = [inds_joined; tpm_sbn_ind];
plim = [min(inds_joined(:)) max(inds_joined(:))]; % Probability colourbar

figure;

% Empirical TPM
pos = [1 3 5 7];
for state = 1 : size(tpm_sbn, 1)
    subplot(size(tpm_sbn, 1), 2, pos(state));
    bar(tpm_sbn(state, :));
    ylim([0.42 0.58]);
end

% Independent TPM
pos = [2 4 6 8];
for state = 1 : size(tpm_sbn_ind, 1)
    subplot(size(tpm_sbn_ind, 1), 2, pos(state));
    bar(tpm_sbn_ind(state, :));
    ylim([0.42 0.58]);
end


%% Composition
% Part d) of the figure

% Get composition
% Order of concepts - 0, 1, 01
unpart = [...
    phi.big_mips{1}.unpartitioned_constellation{1}.phi...
    phi.big_mips{1}.unpartitioned_constellation{2}.phi...
    phi.big_mips{1}.unpartitioned_constellation{3}.phi...
    ];
part = [...
    phi.big_mips{1}.partitioned_constellation{1}.phi...
    phi.big_mips{1}.partitioned_constellation{2}.phi...
    phi.big_mips{1}.partitioned_constellation{3}.phi...
    ];
diff = unpart - part;

% Set of line-sources (each corresponds to a concept)
% Each contains a list of destination concepts
sources{1} = [3];
sources{2} = [3];
sources{3} = [];

concept_labels = {'A', 'B', 'AB'};

% Colour per concept-order
colours = [1 1 2];

% x-axis is concept order
x = [1 1 2];
% y-axis is arbitrary
y = [1 2 1.5]; % evenly space 2nd order concept between first orders
[ysorted, yidx] = sort(y);
% z-axis is small-phi

fig = figure;
set(gcf, 'Position', [0 0 2100/2 250]);
set(gcf, 'Color', [1 1 1]);
set(gcf, 'InvertHardCopy', 'off'); % For keeping the black background when printing
set(gcf, 'RendererMode', 'manual');
set(gcf, 'Renderer', 'painters');

colormap prism

size = 100;
margin = 0.2;

% Establish axes
subplot(1, 4, 1);
scatter3(x, y, unpart, size, colours', 'o', 'filled'); % actual values
hold on;
for source = 1:length(sources)
    for dest = sources{source}
        % Draw lines
        line([x(source) x(dest)], [y(source) y(dest)], [unpart(source) unpart(dest)], 'Color', 'k');
    end
end
scatter3(x, y, [0 0 0], size, colours', 'o'); % 'shadow' values
for source = 1:length(sources)
    for dest = sources{source}
        % Draw lines
        line([x(source) x(dest)], [y(source) y(dest)], 'Color', 'k', 'Linestyle', ':');
    end
end
xlim([min(x)-margin max(x)+margin]); ylim([min(y)-margin max(y)+margin]); zlim([0.003 0.1]);
view([135 30]);
zlabel('\phi', 'rotation', 0);
set(gca, 'ZTick', [0.01 0.07], 'YTick', [1 2], 'YTickLabel', {'A', 'B'}, 'XTick', 2, 'XTickLabel', 'AB');
set(gca, 'zscale', 'log');
box on
axis square

% Unpartitioned comp
subplot(1, 4, 2);
scatter3(x, y, unpart, size, colours', 'o', 'filled'); set(gca, 'zscale', 'log');
xlim([min(x)-margin max(x)+margin]); ylim([min(y)-margin max(y)+margin]); zlim([0.003 0.1]);
view([90 0]);
% Draw lines
for source = 1:length(sources)
    for dest = sources{source}
        line([x(source) x(dest)], [y(source) y(dest)], [unpart(source) unpart(dest)], 'Color', 'k');
    end
end
grid off;
title('UP');
set(gca, 'YTick', ysorted, 'YTickLabel', concept_labels(yidx));
set(gca, 'ZTick', sort(unpart));
box on
axis square

% Partitioned comp
subplot(1, 4, 3);
scatter3(x, y, part, size, colours', 'o', 'filled');
xlim([min(x)-margin max(x)+margin]); ylim([min(y)-margin max(y)+margin]); zlim([0.003 0.1]);
view([90 0]);
% Draw lines
for source = 1:length(sources)
    for dest = sources{source}
        line([x(source) x(dest)], [y(source) y(dest)], [part(source) part(dest)], 'Color', 'k');
    end
end
hold on; scatter3(x(3), y(3), unpart(3), size, 'o', 'filled');
grid off;
title('P');
set(gca, 'YTick', ysorted, 'YTickLabel', concept_labels(yidx));
set(gca, 'ZTick', sort(part));
set(gca, 'zscale', 'log');
box on
axis square

% Diff comp
subplot(1, 4, 4);
scatter3(x, y, diff, size, colours', 'o', 'filled');
xlim([min(x)-margin max(x)+margin]); ylim([min(y)-margin max(y)+margin]);
view([90 0]);
% Draw lines
for source = 1:length(sources)
    for dest = sources{source}
        line([x(source) x(dest)], [y(source) y(dest)], [diff(source) diff(dest)], 'Color', 'k');
    end
end
grid off;
title('UP-P');
set(gca, 'YTick', ysorted, 'YTickLabel', concept_labels(yidx));
box on
axis square

%% Print

figure_name = 'figures/fig1b';

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG

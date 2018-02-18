%{
Plots aggregrated main effects. Statistical tests of effects are done by
the LME simulated likelihood ratio tests.

Plotted effects:
    nChannels
    tau
    condition
%}

%% Setup

flies = (1:13);
nChannels = (2:4);
taus = [4 8 16];
conditions = (1:2);

ylabel_text = struct();
ylabel_text.three = '\Phi';
ylabel_text.star = '\Phi*';

condition_shapes = 'ox';
condition_offset = [-0.1 0.1];
tau_colours = 'rgb';

%% Load

phi_type = 'three'; % 'three' or 'star'

data_nChannels = '2t4';
data_detrended = 0;
data_zscored = 0;

data_directory = 'results/';
data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim'...
    '_detrend' num2str(data_detrended)...
    '_zscore' num2str(data_zscored)...
    '_nChannels' data_nChannels...
    '_phi' phi_type...
    ];

results_directory = 'analysis_results/';
results_filename = [data_filename '_lmeStats'];

disp('loading');

load([results_directory results_filename '.mat']);

disp('loaded');

%% Un-log phi values

phi_values = phi_table(:, 7);
phi_values = table2array(phi_values);
phi_values = exp(phi_values);
phi_values = array2table(phi_values);
phi_table(:, 7) = phi_values;

%% Get data to plot



% Average and std across flies

plot_data = struct();
plot_data.means = zeros(length(conditions), length(taus), length(nChannels));
plot_data.stds = zeros(length(conditions), length(taus), length(nChannels));

for nChannels_counter = 1 : length(nChannels)
    for tau_counter = 1 : length(taus)
        for condition_counter = 1 : length(conditions)
            rows = phi_table(...
                phi_table.nChannels == nChannels(nChannels_counter) &...
                phi_table.tau == tau_counter &...
                phi_table.condition == conditions(condition_counter),...
                :);
            % Average sets within flies
            fly_data = zeros(length(flies), 1);
            for fly_counter = 1 : length(flies)
                fly_rows = rows(rows.fly == flies(fly_counter), 7);
                fly_data(fly_counter) = mean(table2array(fly_rows));
            end
            
            % Average across flies
            plot_data.means(condition_counter, tau_counter, nChannels_counter) = mean(fly_data);
            plot_data.stds(condition_counter, tau_counter, nChannels_counter) = std(fly_data) / sqrt(length(flies));
        end
    end
end

%% Plot
% Scheme: one subplot per condition
% Within each condition, plot for each nChannels and tau

figure;

nChannels_xpos = [2 5 8];
nChannels_ypos = -0.005;
figure_compress = 0.05;

x_pos = 1;
for nChannels_counter = 1 : length(nChannels)
    for tau_counter = 1 : length(taus)
        for condition_counter = 1 : length(conditions)
            errorbar(x_pos + condition_offset(condition_counter),...
                plot_data.means(condition_counter, tau_counter, nChannels_counter),...
                plot_data.stds(condition_counter, tau_counter, nChannels_counter),...
                ['k' condition_shapes(condition_counter)]);
            hold on;
            axis([0 10 0 0.03]);
        end
        x_pos = x_pos + 1;
    end
    % nChannels grouping xlabels
    h = text(nChannels_xpos(nChannels_counter), nChannels_ypos, [num2str(nChannels(nChannels_counter)) 'ch'], 'Fontweight', 'bold');
    pos = get(h, 'Position');
    ext = get(h, 'Extent');
    pos(1) = pos(1) - ext(3)/2;
    set(h, 'Position', pos);
end

pos = get(gca, 'position');
pos(2) = pos(2) + figure_compress;
set(gca, 'Position', pos);

set(gca,...
    'XTick', (1:length(taus)*length(nChannels)),...
    'XTickLabel', [4 8 16 4 8 16 4 8 16],...
    'YTick', [0 0.015 0.03]);

%% Plot
% Scheme: treat nChannels as independent, so show how each condition
% affects each tau condition within each nChannels

figure;

for nChannels_counter = 1 : length(nChannels)
    subplot(1, length(nChannels), nChannels_counter);
    bar(plot_data.means(:, :, nChannels_counter));
    axis([0.5 2.5 0 0.025]);
    title([num2str(nChannels(nChannels_counter)) ' channels']);
end

% Plot for 1 nChannels, (to show trend)


%% Plot
% Scheme: treat taus as indepedent, so show how each condition affects each
% nChannels condition within each tau

plot_data_perTau = struct();
plot_data_perTau.means = permute(plot_data.means, [1 3 2]);
plot_data_perTau.stds = permute(plot_data.stds, [1 3 2]);

figure;

for tau_counter = 1 : length(taus)
    s = subplot(1, length(taus), tau_counter);
    [bars, errors] = barwitherr(plot_data_perTau.stds(:, :, tau_counter), plot_data_perTau.means(:, :, tau_counter));
    axis([0.5 2.5 0 0.028]);
    title(['\tau = ' num2str(taus(tau_counter)) ' ms'], 'FontSize', 15);
    ylabel(char(981), 'FontSize', 15, 'rotation', 0);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'air', 'iso'}, 'YTick', [0 0.01 0.02]);
    s.XAxis.FontSize = 15;
    set(errors, 'LineWidth', 0.1, 'CapSize', 0);
    set(bars(3), 'FaceColor', 'g');
    legend('2ch', '3ch', '4ch');
end



%% Focus on main effect of condition

%{
Purpose of this plot is to show that phi decreases for ALL flies

For each fly, size, average across all channel sets

x axis must be 13 flies
y axis should be phi
shape/colour axis can be condition/tau
errorbars can be std across channel sets
subplot axis can be number of channels
%}

plot_tau = 2;
plot_nChannels = 3;

% Obtain data and average across sets
plot_data = struct();
plot_data.means = zeros(length(flies), length(conditions), length(taus), length(nChannels));
plot_data.stds = zeros(length(flies), length(conditions), length(taus), length(nChannels));
for nChannels_counter = 1 : length(nChannels)
    for tau_counter = 1 : length(taus)
        for condition_counter = 1 : length(conditions)
            for fly_counter = 1 : length(flies)
                rows = phi_table(...
                    phi_table.nChannels == nChannels(nChannels_counter) &...
                    phi_table.tau == tau_counter &...
                    phi_table.condition == conditions(condition_counter) &...
                    phi_table.fly == fly_counter,...
                    :);
                % Average across sets
                phi_values = table2array(rows(:, 7));
                %phi_values_iso = table2array(rows_iso(:, 7));
                plot_data.means(fly_counter, condition_counter, tau_counter, nChannels_counter) = mean(phi_values);
                plot_data.stds(fly_counter, condition_counter, tau_counter, nChannels_counter) = std(phi_values);
            end
        end
    end
end

figure;
for nChannels_counter = 1 : length(nChannels)
    subplot(1, length(nChannels), nChannels_counter);
    for tau_counter = 1:length(taus)
        for fly_counter = 1 : length(flies)
            for condition_counter = 1 : length(conditions)
%                 errorbar(fly_counter + condition_offset(condition_counter),...
%                     plot_data.means(fly_counter, condition_counter, tau_counter, nChannels_counter),...
%                     plot_data.stds(fly_counter, condition_counter, tau_counter, nChannels_counter),...
%                     ['k' condition_shapes(condition_counter)]); hold on;
                scatter(fly_counter,...
                    plot_data.means(fly_counter, condition_counter, tau_counter, nChannels_counter),...
                    ['k' condition_shapes(condition_counter)]); hold on;
                axis([0 14 0 0.08]);
                set(gca, 'XTick', (1:13), 'XTickLabel', [], 'YTick', [0 0.04 0.08]);
                xlabel('fly'); ylabel(ylabel_text.(phi_type));
            end
        end
    end
end

% Stats test
ps = zeros(length(nChannels), length(taus));
stats = cell(length(nChannels), length(taus));
for nChannels_counter = 1 : length(nChannels)
    for tau_counter = 1 : length(taus)
        air = plot_data.means(:, 1, tau_counter, nChannels_counter);
        iso = plot_data.means(:, 2, tau_counter, nChannels_counter);
        [ps(nChannels_counter, tau_counter), ~, stats{nChannels_counter, tau_counter}]= signrank(air, iso);
    end
end

% Bonferroni correction
ps = ps * numel(ps);

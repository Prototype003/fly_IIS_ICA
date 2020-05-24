%% Description

%{
Plots big-phi at each tau level
%}

%% Setup

load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICAAllTrials_nComponents4_phithree_nChannels4_globalTPM0.mat');

%% Plot big-phi
% 1 plot per fly, comparing all trials of wake vs anesthesia

cond_colours = {'r', 'b'};

big_phis = phis{1}.phis;

figure;
for fly = 1 : size(big_phis, 3)
    subplot(4, 4, fly);
    for condition = 1 : size(big_phis, 4)
        plot(phis{1}.taus, squeeze(big_phis(:, :, fly, condition, :)), cond_colours{condition});
        hold on;
    end
    set(gca, 'XScale', 'log');
    
    xlabel('tau-step (ms)');
    ylabel('\Phi');
    title(['fly' num2str(fly)]);
    
    ylim([min(big_phis(:)) max(big_phis(:))]);
end
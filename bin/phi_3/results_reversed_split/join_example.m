trials = 8;
conds = 2;

vals_reversed = zeros(trials, conds);

for cond = 1 : conds
    for trial = 1 : trials
        filename =...
            ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStimreversed'...
            '_nChannels4_globalTPM0_'...
            'f01c' num2str(cond) 'tau4s0876t' num2str(trial) '.mat'...
            ];
        tmp = load(filename);
        vals_reversed(trial, cond) = tmp.phi.phi;
    end
end

%% original time order

load('../results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM0.mat');

channel_set = 876;
fly = 1;

phi_vals = phis{1}.phis;

vals_orig = squeeze(phi_vals(channel_set, :, fly, :));

%% Plot

figure;

plot(vals_orig(:, 1), 'r-'); hold on;
plot(vals_orig(:, 2), 'b-');
plot(vals_reversed(:, 1), 'r--');
plot(vals_reversed(:, 2), 'b--');

ylabel('\Phi'); xlabel('trial');

legend('wake', 'anest', 'w\_rev', 'a\_rev', 'location', 'northwest');
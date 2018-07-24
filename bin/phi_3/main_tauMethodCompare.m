%% Description

%{
Compares phi values for binned tau values

for fly 1:
    2ch : 61
    4ch : 1036
%}

%% Setup

nChannels = 4;
channel_set = 1036;


taus = [1 4 8 12 16 24 32 48 64];
%taus = (1:32);
source_dir = 'results_split/';

source_file_prefix = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels' num2str(nChannels) '_globalTPM1'];
source_file_suffix = 't1.mat';

%% 

tau_mat = zeros(16, 16);

taus = [1 4 8 16];
tau_string = 'tau';
for tau = taus
    infix = ['_f01c1' tau_string num2str(tau) 's' num2str(channel_set) source_file_suffix];
    load([source_dir source_file_prefix infix]);
    tau_mat(tau, 1) = phi.phi;
end



taus = [2 3 4 8 16];
tau_string = 'tauBin';
for tau = taus
    infix = ['_f01c1' tau_string num2str(tau) 's' num2str(channel_set) source_file_suffix];
    load([source_dir source_file_prefix infix]);
    tau_mat(tau, tau) = phi.phi;
end

%%

figure;

imagesc(tau_mat);

xlabel('bin size (ms)');
ylabel('real tau (after binning) (ms)');

c = colorbar;
ylabel(c, '\phi');

title('fly 1, wake, 4ch set 1036');
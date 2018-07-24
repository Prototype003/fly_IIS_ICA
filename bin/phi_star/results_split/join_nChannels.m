%% Description

%{

Joins each nChannels result

%}

%% Setup

phi_type = 'SI'; % 'star' or 'SI';
global_type = 'nonGlobal'; % 'global' or 'nonGlobal'

nChannels = (2:4);
nChannels_string = [num2str(nChannels(1)) 't' num2str(nChannels(end))];

source_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels';
source_suffix = ['_phi' phi_type '_' global_type '.mat'];

out_dir = '../results/';
out_file = [source_prefix nChannels_string source_suffix];

%% Join

phis = cell(0);
for nChannels_counter = 1 : length(nChannels)
    source_file = [source_prefix num2str(nChannels(nChannels_counter)) source_suffix];
    
    tmp = load(source_file);
    
    phis{nChannels_counter} = tmp.phis{1};
    phis{nChannels_counter}.channel_sets = nchoosek((1:15), nChannels(nChannels_counter));
end

%% Save

save([out_dir out_file], 'phis');
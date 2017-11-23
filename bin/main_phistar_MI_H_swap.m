%% Description

%{

WARNING: OVERWRITES DATAFILES

For phi-star results:

Swap struct fields mis with hs, and mi_stars with hconds

%}

%% Setup

flies = (1:13);

% Add these variables
nChannels_vars = [2 3 4];
taus = [4 8 16];
nBins = 1;

data_directory = 'results/preformatted_results/';
prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_medianSplit0_phistar_allPartitionsfly';
suffix = '.mat';

%% Swap

for fly = 1 : length(flies)
    disp(fly);
    
    % Load
    phis = load([data_directory prefix num2str(fly) suffix]);
    phis = phis.phis;
    
    for nChannels = 1 : length(phis)
        % Swap MIs with Hs
        mis = phis{nChannels}.partitions_hs;
        mi_stars = phis{nChannels}.partitions_hconds;
        
        phis{nChannels}.partitions_hs = phis{nChannels}.partitions_mis;
        phis{nChannels}.partitions_hconds = phis{nChannels}.partitions_mi_stars;
        
        phis{nChannels}.partitions_mis = mis;
        phis{nChannels}.partitions_mi_stars = mi_stars;
        
        phis{nChannels}.nChannels = nChannels_vars(nChannels);
        phis{nChannels}.taus = taus;
        phis{nChannels}.nBins = nBins;
    end
    
    % Save
    save([data_directory prefix num2str(fly) suffix], 'phis');
    
end
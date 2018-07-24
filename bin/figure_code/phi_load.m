function [phis, measure_string] = phi_load(measure, global_tpm, root_directory)
%PHI_LOAD
% Loads phi results (either for phi-3 or phi-star)
%
% Inputs:
%   measure: string
%       'phi_three' or 'phi_star' or 'phi_SI'
%   global_tpm: integer
%       0 - get phi values using 2.25s TPM/covariance observation (per trial)
%       1 - get phi values computed 18s TPM/covariance observation (across trials)
%   root_directory: string
%       Specifies location of the root folder (fly_phi/bin/) of the project
%       Used for locating results directory locations
%
% Outputs:
%   phis: cell vector
%       Each cell contains phi results for one network size
%   measure_string: string
%       A string which can be used as a label when making figures

if strcmp(measure, 'phi_three') % Load phi-three results
    measure_string = '\Phi^{3.0}';
    
    if global_tpm == 1 % Global TPM
        disp('loading');
        data_directory = [root_directory 'results/'];
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree_global.mat'];
        load([data_directory data_filename]);
        
        % Fix python indexing and rename measure specific variable names to generic names
        for nChannels_counter = 1 : length(phis)
            phis{nChannels_counter}.channel_sets = double(phis{nChannels_counter}.channel_sets) + 1;
            phis{nChannels_counter}.phis = phis{nChannels_counter}.phi_threes;
            rmfield(phis{nChannels_counter}, 'phi_threes');
        end
        
        disp('loaded');
    else % global_tpm == 0 % TPMs per trial
        phis = cell(1, 3);
        
        data_directory = [root_directory 'results/'];
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t2_phithree_nonGlobal.mat'];
        disp('loading 2ch');
        tmp = load([data_directory data_filename]);
        phis{1} = tmp.phis{1};
        
        data_directory = [root_directory 'results/'];
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal.mat'];
        disp('loading 3ch');
        tmp = load([data_directory data_filename]);
        phis{2} = tmp.phis{1};
        
        data_directory = [root_directory 'results/'];
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal.mat'];
        disp('loading 4ch');
        tmp = load([data_directory data_filename]);
        phis{3} = tmp.phis{1};
        clear tmp % Clear up memory
        
        % Fix python indexing and rename measure specific variable names to generic names
        for nChannels_counter = 1 : length(phis)
            if nChannels_counter == 1
                phis{nChannels_counter}.channel_sets = double(phis{nChannels_counter}.channel_sets) + 1;
            end
            phis{nChannels_counter}.phis = phis{nChannels_counter}.phi_threes;
            phis{nChannels_counter} = rmfield(phis{nChannels_counter}, 'phi_threes');
        end
        
%         data_directory = [root_directory 'results/'];
%         data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal_tau4.mat'];
%         disp('loading 3ch');
%         tmp = load([data_directory data_filename]);
%         phis{2} = tmp.phis{1};
%         
%         data_directory = [root_directory 'results/'];
%         data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal_tau4.mat'];
%         disp('loading 4ch');
%         tmp = load([data_directory data_filename]);
%         phis{3} = tmp;
%         clear tmp % Clear up memory
%         
%         % Fix python indexing and rename measure specific variable names to generic names
%         for nChannels_counter = 1 : length(phis)
%             if nChannels_counter == 3
%                 phis{nChannels_counter}.channel_sets = double(phis{nChannels_counter}.phis{1}.channel_sets);
%             else
%                 phis{nChannels_counter}.channel_sets = double(phis{nChannels_counter}.channel_sets) + 1;
%             end
%             phis{nChannels_counter}.phis = phis{nChannels_counter}.phi_threes;
%             phis{nChannels_counter} = rmfield(phis{nChannels_counter}, 'phi_threes');
%         end
        
        disp('loaded');
    end
    
elseif strcmp(measure, 'phi_star') % Load phi-star results
    measure_string = '\Phi*';
%     % Results from calculation using Gaussian assumption
%     if global_tpm == 1 % Global covariance
%         disp('loading');
%         
%         data_directory = [root_directory 'results/'];
%         data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_medianSplit0_phistar_global.mat'];
%         load([data_directory data_filename]);
%         
%         % Rename measure specific variable names to generic names
%         for nChannels_counter = 1 : length(phis)
%             phis{nChannels_counter}.phis = phis{nChannels_counter}.phi_stars;
%             phis{nChannels_counter} = rmfield(phis{nChannels_counter}, 'phi_stars');
%         end
%         
%         disp('loaded');
%         
%     else % global_tpm == 0 % Covariance per trial
%         disp('loading');
%         
%         data_directory = [root_directory 'results/'];
%         data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar.mat'];
%         load([data_directory data_filename]);
%         
%         % Rename measure specific variable names to generic names
%         for nChannels_counter = 1 : length(phis)
%             phis{nChannels_counter}.phis = phis{nChannels_counter}.phi_stars;
%             phis{nChannels_counter} = rmfield(phis{nChannels_counter}, 'phi_stars');
%         end
%         
%         disp('loaded');
%     end
    
    disp('loading');
    phis = cell(3, 1);
    if global_tpm == 1 % 1 big trial
        data_directory = [root_directory 'phi_star/results/'];
        data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar_global.mat';
        tmp = load([data_directory data_filename]);
        phis = tmp.phis;
    else % global_tpm == 0 % 8 trials
        data_directory = [root_directory 'phi_star/results/'];
        data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar_nonGlobal.mat';
        tmp = load([data_directory data_filename]);
        phis = tmp.phis;
%         data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2_phistar_nonGlobal.mat';
%         tmp = load([data_directory data_filename]);
%         phis{1} = tmp.phis{1};
%         phis{1}.channel_sets = nchoosek((1:15), 2);
%         data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3_phistar_nonGlobal.mat';
%         tmp = load([data_directory data_filename]);
%         phis{2} = tmp.phis{1};
%         phis{2}.channel_sets = nchoosek((1:15), 3);
%         data_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4_phistar_nonGlobal.mat';
%         tmp = load([data_directory data_filename]);
%         phis{3} = tmp.phis{1};
%         phis{3}.channel_sets = nchoosek((1:15), 4);
    end
    disp('loaded');
    
elseif strcmp(measure, 'phi_SI')
    disp('loading');
    
    measure_string = '\Phi^{SI}';
    
    if global_tpm == 1
        data_directory = [root_directory 'phi_star/results/'];
        data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phiSI_global.mat';
        
        tmp = load([data_directory data_file]);
        phis = tmp.phis;
    else % global_tpm == 0 % 8 trials
        data_directory = [root_directory 'phi_star/results/'];
        data_file = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phiSI_nonGlobal.mat';
        
        tmp = load([data_directory data_file]);
        phis = tmp.phis;
    end
    
    disp('loaded');
    
end


end


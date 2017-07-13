%% DESCRIPTION

%{
This script joins the per-fly results from phi3 calculation -
    Calculations were split into nChannels (2:3) and (4)
    For nChannels 4, calculations were split by fly
%}

%% SETUP

prep_detrend = 0;
prep_zscore = 0;

flies = (1:13) - 1; % Because of python 0-based indexing

nFlies = length(flies);
nChannels_base = (2:3);
nChannels_subsequent = (4);
nChannels = [nChannels_base nChannels_subsequent];

data_directory = 'results/preformatted_results/';
data_filename_common = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim' '_detrend' num2str(prep_detrend) '_zscore' num2str(prep_zscore)];

results_directory = 'results/';
results_filename = [data_filename_common '_nChannels' num2str(nChannels_base(1)) 't' num2str(nChannels_subsequent(end)) '_phithree'];

%% Load base results

nChannels_string = ['_nChannels' num2str(nChannels_base(1)) 't' num2str(nChannels_base(end))];
result_tmpChannel = load([data_directory data_filename_common nChannels_string '_phithree' '.mat']);

phis = result_tmpChannel.phis;

%% Reformat MIP for base results

% Order of appending in python script was:
% for fly, for condition, for tau, for channel-set, for state-index

for nChannels_counter = 1 : length(phis)
    channels_used = phis{nChannels_counter}.nChannels;
    
    % 4th dimension is nFlies (should be the same as specified during setup)
    [nStates, nSets, ~, nConditions, nTaus] = size(phis{nChannels_counter}.state_phis);
    
    % Looks like the partitions only splits into 2 groups
    % MIPs are calculated per state, not per trial
    % (only variance across trials is the state-weighting when averaging phi)
    mips_formatted = cell(nStates, nSets, nFlies, nConditions, nTaus);
    
    % This is for the case where nChannels = 1 (where partitions always consist of 1 element and thus MIPs were stored in a matrix instead of cell)
    if isa(phis{nChannels_counter}.mips, 'int64')
        phis{nChannels_counter}.mips = num2cell(phis{nChannels_counter}.mips);
    end
    
    % Go through each fly, condition, tau, channel set, and state in the same order as they were iterated through in phi-3 python script
    % Go through mips in a linear fashion (each iteration's mip was appended to a 1D list)
    mip_counter = 1;
    for fly = 1 : nFlies
        for condition = 1 : nConditions
            for tau = 1 : nTaus
                for set = 1 : nSets
                    for state = 1 : nStates
                        cut{1} = phis{nChannels_counter}.mips(mip_counter, :);
                        mips_formatted(state, set, fly, condition, tau) = cut;
                        mip_counter = mip_counter + 1;
                    end
                end
            end
        end
    end
    
    % Replace the unformatted .mips with the new format
    phis{nChannels_counter}.mips = mips_formatted;
    
end

%% Append other channels and flies

for nChannel_counter = 1 : length(nChannels_subsequent)
    nChannels_string = ['_nChannels' num2str(nChannels_subsequent(nChannel_counter)) 't' num2str(nChannels_subsequent(nChannel_counter))];
    
    % Load first fly to append subsequent flies to
    fly_counter = 1;
    fly = flies(fly_counter);
    disp(fly);
    result_tmpChannel = load([data_directory data_filename_common nChannels_string '_phithreefly' num2str(fly) '.mat']);
    result_tmpChannel = result_tmpChannel.phis;
    for channels_counter = 1 : length(result_tmpChannel)
        % Reformat mips into cell array with appropriate dimensions for the first fly
        result_tmpChannel{channels_counter}.mips = format_mips(result_tmpChannel{channels_counter}.mips, size(result_tmpChannel{channels_counter}.state_phis));
        
        % Load and append subsequent flies
        for fly_counter = 2 : nFlies
            fly = flies(fly_counter);
            disp(fly);
            result_tmpFly = load([data_directory data_filename_common nChannels_string '_phithreefly' num2str(fly) '.mat']);
            result_tmpFly = result_tmpFly.phis;
            
            % Reformat mips into cell array
            result_tmpFly{channels_counter}.mips = format_mips(result_tmpFly{channels_counter}.mips, size(result_tmpFly{channels_counter}.state_phis));
            
            % Concatenate variables
            result_tmpChannel{channels_counter}.phi_threes = cat(3, result_tmpChannel{channels_counter}.phi_threes, result_tmpFly{channels_counter}.phi_threes);
            result_tmpChannel{channels_counter}.state_counters = cat(4, result_tmpChannel{channels_counter}.state_counters, result_tmpFly{channels_counter}.state_counters);
            result_tmpChannel{channels_counter}.state_phis = cat(3, result_tmpChannel{channels_counter}.state_phis, result_tmpFly{channels_counter}.state_phis);
            result_tmpChannel{channels_counter}.tpms = cat(4, result_tmpChannel{channels_counter}.tpms, result_tmpFly{channels_counter}.tpms);
            result_tmpChannel{channels_counter}.mips = cat(3, result_tmpChannel{channels_counter}.mips, result_tmpFly{channels_counter}.mips);
        end
        
        % Append fly-joined channel-set results
        phis{length(phis)+1} = result_tmpChannel{channels_counter};
    end
    
end

%% SAVE
disp('Saving');
if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename '.mat'], 'phis');

%% Reformatting of .mips structure

function [mips_formatted] = format_mips(mips, state_phi_dimensions)
%
% Inputs:
%   mips = is a 2D cell array (or 2D matrix) for a single fly - each row is a MIP,
%       unidirectional connections from column 1 to column 2 are those which are cut in the MIP
%   state_phi_dimensions = vector, output of size(.state_phis)
%
% Outputs:
%   mips_formatted = cell array with dimensions as specified by state_phi_dimensions, with an additional dimension
%       (for holding the cut information)

dimensions_cell = num2cell(state_phi_dimensions);
[nStates, nSets, nFlies, nConditions, nTaus] = dimensions_cell{:};

mips_formatted = cell([nStates, nSets, nFlies, nConditions, nTaus]);

mip_counter = 1;
for fly = 1 : nFlies
    for condition = 1 : nConditions
        for tau = 1 : nTaus
            for set = 1 : nSets
                for state = 1 : nStates
                    cut{1} = mips(mip_counter, :);
                    mips_formatted(state, set, fly, condition, tau) = cut;
                    mip_counter = mip_counter + 1;
                end
            end
        end
    end
end

end
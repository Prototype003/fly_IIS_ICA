%% DESCRIPTION

%{
This script joins the per-fly results from phi3 calculation for all partitions -
    Calculations were split into nChannels (2), (3), and (4)
    For nChannels 4, calculations were split by fly

Note:
Partitions were saved in the following manner:
    If number of channels between groups is unequal: cell array,
    partition{1} --/--> partition{2}
    If number of channels between groups is equal: int64 array,
    partition(1, :) --/--> partition(2, :)
%}

%% SETUP

prep_detrend = 0;
prep_zscore = 0;

flies = (1:13) - 1; % Because of python 0-based indexing

nFlies = length(flies);
nChannels_base = (2);
nChannels_subsequent = (3);
nChannels_subsequent_flySplitted = (4);
nChannels = [nChannels_base nChannels_subsequent nChannels_subsequent_flySplitted];

suffix_string = '_phithree_allPartitions';

data_directory = 'results/preformatted_results/';
data_filename_common = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim' '_detrend' num2str(prep_detrend) '_zscore' num2str(prep_zscore)];

results_directory = 'results/';
results_filename = [data_filename_common '_nChannels' num2str(nChannels_base(1)) 't' num2str(nChannels(end)) suffix_string];

%% Load base results

nChannels_string = ['_nChannels' num2str(nChannels_base(1)) 't' num2str(nChannels_base(end))];
result_tmpChannel = load([data_directory data_filename_common nChannels_string suffix_string '.mat']);

phis = result_tmpChannel.phis;

%% Reformat partitions for base results

for nChannels_counter = 1 : length(phis)
    
%     for mip_counter = 1 : numel(phis{nChannels_counter}.mips)
%         if isa(phis{nChannels_counter}.mips{mip_counter}, 'int64')
%             phis{nChannels_counter}.mips{mip_counter} = format_partition(phis{nChannels_counter}.mips{mip_counter});
%         end
%     end
%     
%     for partition_counter = 1 : numel(phis{nChannels_counter}.state_partitions)
%         if isa(phis{nChannels_counter}.state_partitions{partition_counter}, 'int64')
%             phis{nChannels_counter}.state_partitions{partition_counter} = format_partition(phis{nChannels_counter}.state_partitions{partition_counter});
%         end
%     end
    
end

%% Append other nChannels which are not split by flies

for nChannel_counters = 1 : length(nChannels_subsequent)
    nChannels_string = ['_nChannels' num2str(nChannels_subsequent(nChannel_counters)) 't' num2str(nChannels_subsequent(nChannel_counters))];
    disp(['nChannels ' num2str(nChannels_subsequent(nChannel_counters))]);
    
    result_tmpChannel = load([data_directory data_filename_common nChannels_string suffix_string '.mat']);
    result_tmpChannel = result_tmpChannel.phis;
    
%     for mip_counter = 1 : numel(result_tmpChannel{nChannels_counter}.mips)
%         if isa(result_tmpChannel{nChannels_counter}.mips{mip_counter}, 'int64')
%             result_tmpChannel{nChannels_counter}.mips{mip_counter} = format_partition(result_tmpChannel{nChannels_counter}.mips{mip_counter});
%         end
%     end
%     
%     for partition_counter = 1 : numel(result_tmpChannel{nChannels_counter}.state_partitions)
%         if isa(result_tmpChannel{nChannels_counter}.state_partitions{partition_counter}, 'int64')
%             result_tmpChannel{nChannels_counter}.state_partitions{partition_counter} = format_partition(result_tmpChannel{nChannels_counter}.state_partitions{partition_counter});
%         end
%     end
    
    phis{length(phis)+1} = result_tmpChannel{nChannels_counter};
    
end

%% Append other nChannels which are split by flies

for nChannel_counters = 1 : length(nChannels_subsequent_flySplitted)
    nChannels_string = ['_nChannels' num2str(nChannels_subsequent_flySplitted(nChannel_counters)) 't' num2str(nChannels_subsequent_flySplitted(nChannel_counters))];
    disp(['nChannels ' num2str(nChannels_subsequent_flySplitted(nChannel_counters))]);
    
    % First fly
    fly_counter = 1;
    fly = flies(fly_counter);
    disp(['fly' num2str(fly)]);
    
    % Load
    result_tmpChannel = load([data_directory data_filename_common nChannels_string suffix_string 'fly' num2str(fly) '.mat']);
    result_tmpChannel = result_tmpChannel.phis;
    
%     % Reformat partitions to cell arrays
%     for mip_counter = 1 : numel(result_tmpChannel{nChannels_counter}.mips)
%         if isa(result_tmpChannel{nChannels_counter}.mips{mip_counter}, 'int64')
%             result_tmpChannel{nChannels_counter}.mips{mip_counter} = format_partition(result_tmpChannel{nChannels_counter}.mips{mip_counter});
%         end
%     end
%     for partition_counter = 1 : numel(result_tmpChannel{nChannels_counter}.state_partitions)
%         if isa(result_tmpChannel{nChannels_counter}.state_partitions{partition_counter}, 'int64')
%             result_tmpChannel{nChannels_counter}.state_partitions{partition_counter} = format_partition(result_tmpChannel{nChannels_counter}.state_partitions{partition_counter});
%         end
%     end
    
    % Subsequent flies
    for fly_counter = 2 : nFlies
        fly = flies(fly_counter);
        disp(['fly' num2str(fly)]);
        
        % Load
        result_tmpFly = load([data_directory data_filename_common nChannels_string suffix_string 'fly' num2str(fly) '.mat']);
        result_tmpFly = result_tmpFly.phis;
        
%         % Reformat partitions to cell arrays
%         for mip_counter = 1 : numel(result_tmpFly{nChannels_counter}.mips)
%             if isa(result_tmpFly{nChannels_counter}.mips{mip_counter}, 'int64')
%                 result_tmpFly{nChannels_counter}.mips{mip_counter} = format_partition(result_tmpFly{nChannels_counter}.mips{mip_counter});
%             end
%         end
%         for partition_counter = 1 : numel(result_tmpFly{nChannels_counter}.state_partitions)
%             if isa(result_tmpFly{nChannels_counter}.state_partitions{partition_counter}, 'int64')
%                 result_tmpFly{nChannels_counter}.state_partitions{partition_counter} = format_partition(result_tmpFly{nChannels_counter}.state_partitions{partition_counter});
%             end
%         end
        
        % Concatenate results to previous fly/flies
        result_tmpChannel{nChannels_counter}.phi_threes = cat(3, result_tmpChannel{nChannels_counter}.phi_threes, result_tmpFly{nChannels_counter}.phi_threes);
        result_tmpChannel{nChannels_counter}.state_counters = cat(4, result_tmpChannel{nChannels_counter}.state_counters, result_tmpFly{nChannels_counter}.state_counters);
        result_tmpChannel{nChannels_counter}.state_phis = cat(3, result_tmpChannel{nChannels_counter}.state_phis, result_tmpFly{nChannels_counter}.state_phis);
        result_tmpChannel{nChannels_counter}.tpms = cat(4, result_tmpChannel{nChannels_counter}.tpms, result_tmpFly{nChannels_counter}.tpms);
        result_tmpChannel{nChannels_counter}.mips = cat(3, result_tmpChannel{nChannels_counter}.mips, result_tmpFly{nChannels_counter}.mips);
        result_tmpChannel{nChannels_counter}.state_partitions = cat(4, result_tmpChannel{nChannels_counter}.state_partitions, result_tmpFly{nChannels_counter}.state_partitions);
        result_tmpChannel{nChannels_counter}.state_partitions_phis = cat(4, result_tmpChannel{nChannels_counter}.state_partitions_phis, result_tmpFly{nChannels_counter}.state_partitions_phis);
    end
    
    phis{length(phis)+1} = result_tmpChannel{nChannels_counter};
end

%% SAVE
disp('Saving');
if ~isdir(results_directory)
    mkdir(results_directory)
end
save([results_directory results_filename '.mat'], 'phis', '-v7.3');

%% Reformat the partition structure

% function [formatted] = format_partition(partition)
% %
% % Inputs:
% %   partition = a 2D int64 array - top row is the partition cut origin,
% %       bottom row is the partition cut destination
% %
% % Outputs:
% %   formatted = a cell array - formatted{1} is the partition cut origin,
% %       formatted{2} is the partition cut destination
% 
% formatted = cell(1, 2);
% formatted{1} = partition(1, :);
% formatted{2} = partition(2, :);
% 
% end

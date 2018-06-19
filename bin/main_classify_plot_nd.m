%% Description

%{
Plot values as a function of channel location.

Plots are 2D for 2-channels, 3D for 3-channels, and 4D for 4-channels
%}

% First couple of sections are just for demonstrating the concept of mapping to 2D space using PCA and tSNE

% %% Attempt to map N-D simplex space to 2D space using PCA
% 
% channels = (1:15);
% 
% nChannels = 2;
% channel_sets = nchoosek(channels, nChannels);
% 
% [coeff, score, latent, ~, explained] = pca(channel_sets); %pca([channel_sets (1:size(channel_sets, 1))']);
% 
% figure;
% scatter(score(:, 1), score(:, 2), [], (1:size(channel_sets, 1)), '.');
% 
% %% Attempt to map N-D simplex space to 2D space using tsne
% 
% channels = (1:15);
% 
% nChannels = 2;
% nDims = 2;
% channel_sets = nchoosek(channels, nChannels);
% 
% mappedX = tsne(channel_sets, (1:size(channel_sets, 1)), nDims);
% 
% figure;
% scatter(mappedX(:, 1), mappedX(:, 2), [], (1:size(channel_sets, 1)), '.'); colorbar;

%% Settings

measure = 'phi_three'; % 'phi_three' or 'phi_star'
tau = 1; % 1 = 4ms; 2 = 8ms; 3 = 16ms
if tau == 1
    tau_string = '4';
elseif tau == 2
    tau_string = '8';
elseif tau == 3
    tau_string = '16';
end

% 0 for tpms/covariances constructed per trial, and for within fly analyses
% 1 for tpms/covariances constructed across trials, and for across fly analyses
global_tpm = 0;

%% Load phi measurement files

% TODO: Maybe it would be better to stick this whole thing into a separate function somewhere....
if strcmp(measure, 'phi_three') % Load phi-three results
    measure_string = '\Phi3';
    
    if global_tpm == 1 % Global TPM
        disp('loading');
        data_directory = 'results/';
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phithree.mat'];
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
        
        data_directory = 'results/';
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t2_phithree_nonGlobal.mat'];
        disp('loading 2ch');
        tmp = load([data_directory data_filename]);
        phis{1} = tmp.phis{1};
        
        data_directory = 'results/';
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels3t3_phithree_nonGlobal_tau4.mat'];
        disp('loading 3ch');
        tmp = load([data_directory data_filename]);
        phis{2} = tmp.phis{1};
        
        data_directory = 'results_split/';
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels4t4_phithree_nonGlobal_tau4.mat'];
        disp('loading 4ch');
        tmp = load([data_directory data_filename]);
        phis{3} = tmp;
        clear tmp
        
        % Fix python indexing and rename measure specific variable names to generic names
        for nChannels_counter = 1 : length(phis)
            if nChannels_counter == 3
                phis{nChannels_counter}.channel_sets = double(phis{nChannels_counter}.phis{1}.channel_sets);
            else
                phis{nChannels_counter}.channel_sets = double(phis{nChannels_counter}.channel_sets) + 1;
            end
            phis{nChannels_counter}.phis = phis{nChannels_counter}.phi_threes;
            phis{nChannels_counter} = rmfield(phis{nChannels_counter}, 'phi_threes');
        end
        
        disp('loaded');
    end
    
else % strcmp(measure, 'phi_star') % Load phi-star results
    measure_string = '\Phi*';
    
    if global_tpm == 1 % Global covariance
        disp('loading');
        
        data_directory = 'results/';
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_medianSplit0_phistar_global.mat'];
        load([data_directory data_filename]);
        
        % Rename measure specific variable names to generic names
        for nChannels_counter = 1 : length(phis)
            phis{nChannels_counter}.phis = phis{nChannels_counter}.phi_stars;
            phis{nChannels_counter} = rmfield(phis{nChannels_counter}, 'phi_stars');
        end
        
        disp('loaded');
        
    else % global_tpm == 0 % Covariance per trial
        disp('loading');
        
        data_directory = 'results/';
        data_filename = ['split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_detrend0_zscore0_nChannels2t4_phistar.mat'];
        load([data_directory data_filename]);
        
        % Rename measure specific variable names to generic names
        for nChannels_counter = 1 : length(phis)
            phis{nChannels_counter}.phis = phis{nChannels_counter}.phi_stars;
            phis{nChannels_counter} = rmfield(phis{nChannels_counter}, 'phi_stars');
        end
        
        disp('loaded');
    end
    
end

%% Load classification results

results_directory = 'workspace_results/';

if strcmp(measure, 'phi_three')
    if global_tpm == 1
        results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_global_classification_across1.mat';
        load([results_directory results_filename]);
        % Fix python indexing
        for nChannels_counter = 1 : length(accuracies)
            accuracies{nChannels_counter}.channel_sets = accuracies{nChannels_counter}.channel_sets + 1;
        end
    else % global_tpm == 0
        results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phithree_nonGlobal_classification.mat';
        load([results_directory results_filename]);
        % Fix python indexing
        python_indexing = [1 1 0];
        for nChannels_counter = 1 : length(accuracies)
            accuracies{nChannels_counter}.channel_sets = accuracies{nChannels_counter}.channel_sets + python_indexing(nChannels_counter);
        end
    end
else % strcmp(measure, 'phi_star')
    if global_tpm == 1
        results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_global_classification_across1.mat';
        load([results_directory results_filename]);
    else % global_tpm == 0
        results_filename = 'split2250_bipolarRerefType1_lineNoiseRemoved_phistar_nonGlobal_classification.mat';
        load([results_directory results_filename]);
    end
end

%% Get mean values per channel

channels = (1 : max(accuracies{1}.channel_sets(:)));
channel_means = cell(length(phis), 1); % We will save channel means so we can use them for selecting large sets of channels

% Storage of channel means
for nChannels_counter = 1 : length(accuracies)
    channel_sets = accuracies{nChannels_counter}.channel_sets;
    
    % Phi values
    tmp_size = size(phis{nChannels_counter}.phis(:, :, :, :, tau));
    tmp_size(1) = length(channels);
    phis{nChannels_counter}.channel_sums = zeros(tmp_size);
    phis{nChannels_counter}.channel_means = zeros(tmp_size);
    phis{nChannels_counter}.set_counters = zeros(tmp_size);
    
    % Accuracies
    tmp_size = size(accuracies{nChannels_counter}.accuracies);
    tmp_size(1) = length(channels);
    accuracies{nChannels_counter}.channel_sums = zeros(tmp_size);
    accuracies{nChannels_counter}.channel_means = zeros(tmp_size);
    accuracies{nChannels_counter}.set_counters = zeros(tmp_size);
    
    % Sum values for each channel (sum across networks which contain the channel)
    for channel = 1 : max(channel_sets(:))
        for channel_set = 1 : size(channel_sets, 1)
            if any(channel_sets(channel_set, :) == channel)
                
                phis{nChannels_counter}.channel_sums(channel, :, :, :) = phis{nChannels_counter}.channel_sums(channel, :, :, :) + phis{nChannels_counter}.phis(channel_set, :, :, :, tau);
                phis{nChannels_counter}.set_counters(channel, :, :, :) = phis{nChannels_counter}.set_counters(channel, :, :, :) + 1;
                
                accuracies{nChannels_counter}.channel_sums(channel, :) = accuracies{nChannels_counter}.channel_sums(channel, :) + accuracies{nChannels_counter}.accuracies(channel_set, :);
                accuracies{nChannels_counter}.set_counters(channel, :) = accuracies{nChannels_counter}.set_counters(channel, :) + 1;
                
            end
        end
    end
    
    % Average values
    phis{nChannels_counter}.channel_means = phis{nChannels_counter}.channel_sums ./ phis{nChannels_counter}.set_counters;
    accuracies{nChannels_counter}.channel_means = accuracies{nChannels_counter}.channel_sums ./ accuracies{nChannels_counter}.set_counters;
    
    % Store separately for saving
    channel_means{nChannels_counter}.phis = phis{nChannels_counter}.channel_means;
    channel_means{nChannels_counter}.accuracies = accuracies{nChannels_counter}.channel_means;
    
end

% Save channel means
save_location = ['workspace_results/channelMeanPhis_globalTPM' num2str(global_tpm) '_tau' tau_string '.mat'];
save(save_location, 'channel_means');

%% t-SNE and PCA mapping (for consistency across plots)
% We do this separately, first, because t-SNE is an iterative technique
% Because it's iterative, the mapping may change each time (slightly) we use it
% Also, PCA components may flip (because we have a symmetric simplex) each time we run it
% (actually, simplex might not be symmetric across all axes, so this might not be true, but let's just do it here anyway)

% If plotting across measures (i.e. for both phi3 and phi-star), only run this section once to preserve mapping

load_maps = 1; % If 1, then we can load pre-generated maps; if 0, generate and save
view_mapping = 1; % show tSNE maps as they are being generated
map_directory = 'workspace_results/';
map_file = 'simplex_maps.mat';

if load_maps == 1
    load([map_directory map_file]);
else
    disp('Building maps');
    
    % Map storage (nChannels x N mapped dimensions)
    % Each row corresponds a simplex space (row+1 dimensional simplex; row 1 corresponds to 2D simplex)
    % Each column corresponds to the dimensionality of the mapped space (column-D)
    % So (a, b) corresponds to (a+1)-D simplex mapping into b-D
    % e.g. (3, 2) corresponds to 4D simplex mapping into 2D
    reduced_space_max_dimensions = 3; % we can only view up to 3D, but 2D is probably the max we want to plot
    tSNE_map = cell(length(phis), reduced_space_max_dimensions);
    pca_map = cell(size(tSNE_map));
    
    for nChannels_counter = 1 : length(phis)
        channel_sets = phis{nChannels_counter}.channel_sets;
        channels = (1:max(channel_sets(:)));
        
        for reduced_dimension = 1 : size(tSNE_map, 2)
            
            % t-SNE mapping
            if view_mapping == 1 % show maps as they are made
                figure; colormap('jet');
                tSNE_map{nChannels_counter, reduced_dimension} = tsne(channel_sets, (1:size(channel_sets, 1)), reduced_dimension); % Provide labels for plotting
            else % view_mapping == 0 % don't show maps as they are made
                tSNE_map{nChannels_counter, reduced_dimension} = tsne(channel_sets, [], reduced_dimension); % Don't provide labels for plotting
            end
            
            % PCA mapping (these will actually be the same across all "reduced_dimensions" (because no dimension reduction occurs here - the axes are simply rotated)
            [coeff, score, latent, tsquared, explained] = pca(channel_sets);
            pca_map{nChannels_counter, reduced_dimension}.coeff = coeff;
            pca_map{nChannels_counter, reduced_dimension}.score = score;
            pca_map{nChannels_counter, reduced_dimension}.latent = latent;
            pca_map{nChannels_counter, reduced_dimension}.tsquared = tsquared;
            pca_map{nChannels_counter, reduced_dimension}.explained = explained;
            
        end
        
    end
    
    save([map_directory map_file], 'tSNE_map', 'pca_map');
    
    disp('Maps built and saved');
end

% Correlation of path-distance/set-center dimensions to PCs 1 and 2
pc_center_corrs = cell(1, length(accuracies));
pc_distance_corrs = cell(1, length(accuracies));
for nChannels_counter = 1 : length(accuracies)
    
    channel_sets = accuracies{nChannels_counter}.channel_sets;
    
    % Set centers
    centers = mean(channel_sets, 2); % Mean across channels in each set
    [pc_center_corrs{nChannels_counter}.r, pc_center_corrs{nChannels_counter}.p] = corr(centers, pca_map{nChannels_counter, 1}.score(:, 1));
    
    % Set path distances
    distances = channel_set_distances(channel_sets);
    [pc_distance_corrs{nChannels_counter}.r, pc_distance_corrs{nChannels_counter}.p] = corr(distances, pca_map{nChannels_counter, 1}.score(:, 2));
end

%% Plot simplex maps
% TODO
% Plots simplex maps (without any measurement values)

%% Plot settings

flies = (1:13); % Which flies to plot for (only for raw phi values)

measure_type = 'value'; % Raw phi values
%measure_type = 'class'; % classification results

%% Plot raw values on 2D mapping

% nMaps = 2; % PCA and tSNE
% reduced_dimensions = 2;
% 
% figure;
% colormap('jet');
% for nChannels_counter = 1 : length(phis)
%     nChannels = phis{nChannels_counter}.nChannels;
%     
%     if strcmp(measure_type, 'value')
%         
%         values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, flies, :, tau), 2), 3), [1 4 2 3]);
%         
%         plot_values = values(:, 1) - values(:, 2); % Awake - Anest
%         plot_values = values(:, 1)./values(:, 2);
%         %plot_values = values(:, 1); % Awake
%         %plot_values = values(:, 2); % Anest
%         
%         value_sizes = (rescale(plot_values, 0.01, 200));
%         %value_sizes = rescale((plot_values).^5, 0.01, 400); % This changes the scale of values to be non-linear (in order to emphasise big/small values; assumes very small original values)
%         
%     else % strcmp(measure_type, 'class');
%         
%         values = mean(accuracies{nChannels_counter}.accuracies, 2);
%         values(values<58) = NaN;
%         
%         plot_values = values;
%         
%         value_sizes = rescale(plot_values, 0.01, 200);
%         
%     end
%     channel_sets = double(phis{nChannels_counter}.channel_sets);
%     
%     % PCA map
%     subplot(nMaps + 1, length(phis), length(phis)*0 + nChannels_counter);
%     map = pca_map{nChannels_counter, reduced_dimensions}.score;
%     scatter(map(:, 1), map(:, 2), value_sizes, (1:size(channel_sets, 1)), '.');
%     title([num2str(nChannels) 'ch']);
%     
%     % tSNE map
%     subplot(nMaps + 1, length(phis), length(phis)*1 + nChannels_counter);
%     map = tSNE_map{nChannels_counter, reduced_dimensions};
%     scatter(map(:, 1), map(:, 2), value_sizes, (1:size(channel_sets, 1)), '.');
%     
%     % Path distance map
%     subplot(nMaps + 1, length(phis), length(phis)*2 + nChannels_counter);
%     centers = mean(channel_sets, 2); % Set centers - mean across channels in each set
%     distances = channel_set_distances(channel_sets); % Path distances
%     scatter(centers, distances, value_sizes, (1:size(channel_sets, 1)), '.');
%     
% %     % Channel set colorbar;
% %     subplot(nMaps + 1, length(phis), length(phis)*2 + nChannels_counter);
% %     scatter(pca_map{nChannels_counter, reduced_dimensions}.score(:, 1), tSNE_map{nChannels_counter, reduced_dimensions}(:, 1), [], (1:size(channel_sets, 1)), '.'); colorbar;
% %     title('Channel set mapping');
% %     xlabel('PC\_1'); ylabel('tSNE\_1');
% %     axis('tight');
% end

%% Plot raw values on 2D mapping using path distance mapping
% Instead of PCA/tSNE to reduce dimensions, use dimensions: channel set center, and distance among channels in set

figure;
colormap('jet');
nChannels_shapes = '*^+'; sizes = 20;
nChannels_shapes = '...'; sizes = 100;
for nChannels_counter = 1 : length(phis)
    nChannels = phis{nChannels_counter}.nChannels;
    
    channel_sets = phis{nChannels_counter}.channel_sets;
    
    % Set centers
    centers = mean(channel_sets, 2); % Mean across channels in each set
    
    % Set path distances
    distances = channel_set_distances(channel_sets);
    
    % Get values
    if strcmp(measure_type, 'value')
        
        values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, flies, :, tau), 2), 3), [1 4 2 3]);
        
        plot_values = values(:, 1) - values(:, 2); cbar_title = [measure_string '_W - ' measure_string '_A']; % Awake - Anest
        plot_values = values(:, 1)./values(:, 2); cbar_title = [measure_string '_W/' measure_string '_N'];
        plot_values = values(:, 1); cbar_title = [measure_string '_W']; % Awake
        %plot_values = values(:, 2); cbar_title = [measure_string '_N']; % Anest
        
        value_sizes = (rescale(plot_values, 0.01, 200));
        %value_sizes = rescale((plot_values).^5, 0.01, 400); % This changes the scale of values to be non-linear (in order to emphasise big/small values; assumes very small original values)
        
    else % strcmp(measure_type, 'class');
        
        values = mean(accuracies{nChannels_counter}.accuracies, 2); cbar_title = 'class. acc. %';
        %values(values<58) = NaN;
        
        plot_values = values;
        
        value_sizes = rescale(plot_values, 0.01, 200);
        
    end
    
    % Plot
    subplot(1, length(phis), nChannels_counter);
    centers_offset_rand = 0;%(-0.1+(0.1--0.1).*rand(length(centers), 1)); % Values can vary by just under +- half of the minimum difference in centers
    distances_offset_rand = 0;%(-0.4+(0.4--0.4).*rand(length(distances), 1)); % Values can vary by just under +- half of the minimum difference in distances
    scatter(centers+centers_offset_rand, distances+distances_offset_rand, sizes, plot_values, nChannels_shapes(nChannels_counter)); c = colorbar;hold on;
    title([num2str(nChannels) 'Ch']);
    title(c, cbar_title);
    xlabel('set center'); ylabel('set path distance');
    ax = gca; ax.Color = 'k';
    %ylim([0 60]);
    
    % Stats
    [r, p] = corr(distances, plot_values);
    disp([num2str(nChannels_counter+1) 'Ch, distances: r=' num2str(r) ' p=' num2str(p)]);
    [r, p] = corr(centers, plot_values);
    disp([num2str(nChannels_counter+1) 'Ch, centers: r=' num2str(r) ' p=' num2str(p)]);
end
% legend('2ch', '3ch', '4ch');
% ax = gca; ax.Color = 'k';

%% Non-linear colormap
% Source: https://au.mathworks.com/matlabcentral/answers/307318-does-matlab-have-a-nonlinear-colormap-how-do-i-make-one

cMap = jet(256);

dataMax = max(plot_values); % Max of most recent plot values, should be for 4 channel values

dataMin = 0;

centerPoint = 0;

scalingIntensity = 4;


x = 1:length(cMap);

x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);

x = scalingIntensity * x/max(abs(x));


x = sign(x).* exp(abs(x));

x = x - min(x); x = x*511/max(x)+1;

newMap = interp1(x, cMap, 1:512);

%% Plot density-normalised values on center/distance mapping
figure;
set(gcf, 'Position', [0 0 2100/1.5 300]);
colormap('jet');
for nChannels_counter = 1 : length(phis)
    nChannels = phis{nChannels_counter}.nChannels;
    
    channel_sets = phis{nChannels_counter}.channel_sets;
    
    % Get values
    if strcmp(measure_type, 'value')
        values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, flies, :, tau), 2), 3), [1 4 2 3]);
        
        % Select 1
        plot_values = values(:, 1) - values(:, 2); fig_tit = 'diff'; cbar_title = [measure_string '_W - ' measure_string '_A']; % Awake - Anest
        plot_values = values(:, 1)./values(:, 2); fig_tit = 'diff_rel'; cbar_title = [measure_string '_W/' measure_string '_N'];
        %plot_values = values(:, 1); fig_tit = 'wake'; cbar_title = [measure_string '_W']; % Awake
        %plot_values = values(:, 2); fig_tit = 'anest'; cbar_title = [measure_string '_N']; % Anest
        
    else % strcmp(measure_type, 'class');
        values = mean(accuracies{nChannels_counter}.accuracies, 2); fig_tit = 'class'; cbar_title = 'class. acc. %';
        %values(values<58) = NaN;
        plot_values = values;
    end
    
    % Set center / set path distance mapping
    centers = mean(channel_sets, 2); % Mean across channels in each set
    distances = channel_set_distances(channel_sets);
    
    % Find minimum delta among map values
    centers_deltas = pdist2(unique(centers), unique(centers));
    centers_deltas(centers_deltas==0) = NaN; % remove 0 distances
    centers_delta = min(centers_deltas(:));
    distances_deltas = pdist2(unique(distances), unique(distances));
    distances_deltas(distances_deltas==0) = NaN; % remove 0 distances
    distances_delta = min(distances_deltas(:));
    
    % Create space mapping (used for indexing into the mapped space)
    centers_axis = (min(centers) - (2*centers_delta) : centers_delta : max(centers) + (2*centers_delta)); % - and + delta is for padding
    distances_axis = (min(distances) - (2*distances_delta) : distances_delta : max(distances) + (2*distances_delta)); % - and + delta is for padding
    [centers_map, distances_map] = meshgrid(centers_axis, distances_axis);
    
    % Create mapped space
    values_map = zeros(size(centers_map)); % Will sum all values with the same coordinates
    values_map_counter = zeros(size(centers_map)); % Keeps count in each coordinate as to how many values have that coordinate
    
    % Populate mapped space
    for value_counter = 1 : length(plot_values)
        x = find(abs(centers_axis - centers(value_counter)) < 0.00001, 1); % This gives the mapped x location
        y = find(distances_axis == distances(value_counter), 1); % This gives the mapped y location
        values_map(y, x) = values_map(y, x) + plot_values(value_counter); % matrix is (rows, columns), corresponding to (y, x)
        values_map_counter(y, x) = values_map_counter(y, x) + 1;
    end
    
    values_map_plot = values_map ./ values_map_counter; % We will use the NaN values for the black background (and maybe interpolation)
    
    % Interpolate values (linearly - each NaN will turn into the average of the cells immediately adjacent to it, excluding diagonals)
    min_value = min(values_map_plot(:)); % We will remove interpolations which are less than the original min value in the plot (to avoid 'ghost/shadow interpolations');
    values_adjacent = zeros(4, 1);
    for y = 2 : size(values_map_plot, 1) - 1
        for x = 2 : size(values_map_plot, 2) - 1
            if isnan(values_map_plot(y, x))
                values_adjacent(1) = values_map_plot(y, x-1);
                values_adjacent(2) = values_map_plot(y, x+1);
                values_adjacent(3) = values_map_plot(y-1, x);
                values_adjacent(4) = values_map_plot(y+1, x);
                if sum(~isnan(values_adjacent)) > 2
                    %values_adjacent(isnan(values_adjacent)) = 0; % Deal with NaNs
                    %values_adjacent(values_adjacent < min_value) = 0; % Don't include 'shadow' interpolations
                    interpolated = mean(values_adjacent(~isnan(values_adjacent)));
                    if interpolated == 0
                        values_map_plot(y, x) = NaN; % Place NaN instead of 0 to avoid affecting the color scaling
                    else
                        values_map_plot(y, x) = interpolated;
                    end
                end
            end
        end
    end
    %values_map_plot(values_map_plot < min_value) = NaN; % Remove ghost interpolation (which will affect the colorscale)
    
    subplot(1, length(phis), nChannels_counter);
    %imagesc(values_map_plot); c = colorbar;
    plot = pcolor(values_map_plot); c = colorbar;
    set(gca, 'color', [0 0 0]); % black background
    set(plot, 'EdgeColor', 'none'); % remove grid outline
    %axis('xy');
    
    title([num2str(nChannels) 'Ch']);
    title(c, cbar_title);
    xlabel('set center'); ylabel('set path distance');
end

figure_name = ['figures_20180606/' fig_tit];
print(figure_name, '-dpng'); % PNG

%% Plot raw values on 3D mapping
% For 4 channels only (can be done with 3 channels also if you want to)

% flies = (1:13);
% 
% nChannels_counter = 3;
% 
% nChannels = phis{nChannels_counter}.nChannels;
% values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, flies, :, tau), 2), 3), [1 4 2 3]);
% channel_sets = double(phis{nChannels_counter}.channel_sets);
% 
% plot_values = values(:, 1) - values(:, 2); % Awake - Anest
% %plot_values = values(:, 1); % Awake
% %plot_values = values(:, 2); % Anest
% 
% value_sizes = (rescale(plot_values, 0.01, 200));
% value_sizes = rescale((plot_values).^5, 0.01, 400); % This changes the scale of values to be non-linear (in order to emphasise big/small values)
% 
% figure;
% colormap('jet');
% 
% % PCA map
% subplot(1, 2, 1);
% map = pca_map{nChannels_counter, 3}.score;
% scatter3(map(:, 1), map(:, 2), map(:, 3), value_sizes, (1:size(channel_sets, 1)), '.'); colorbar;
% 
% % tSNE map
% subplot(1, 2, 2);
% map = tSNE_map{nChannels_counter, 3};
% scatter3(map(:, 1), map(:, 2), map(:, 3), value_sizes, (1:size(channel_sets, 1)), '.'); colorbar;

%% Plot raw values against 1D space

% flies = (1:13);
% 
% nMaps = 2; % PCA and tSNE
% reduced_dimensions = 1;
% 
% figure;
% colormap('jet');
% for nChannels_counter = 1 : length(phis)
%     nChannels = phis{nChannels_counter}.nChannels;
%     
%     if strcmp(measure_type, 'value')
%         
%         values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, flies, :, tau), 2), 3), [1 4 2 3]);
%         
%         plot_values = values(:, 1) - values(:, 2); % Awake - Anest
%         plot_values = values(:, 1)./values(:, 2);
%         %plot_values = values(:, 1); % Awake
%         %plot_values = values(:, 2); % Anest
%         
%         value_sizes = (rescale(plot_values, 0.01, 200));
%         %value_sizes = rescale((plot_values).^5, 0.01, 400); % This changes the scale of values to be non-linear (in order to emphasise big/small values; assumes very small original values)
%         
%     else % strcmp(measure_type, 'class');
%         
%         values = mean(accuracies{nChannels_counter}.accuracies, 2);
%         values(values<58) = NaN;
%         
%         plot_values = values;
%         
%         value_sizes = rescale(plot_values, 0.01, 200);
%         
%     end
%     channel_sets = double(phis{nChannels_counter}.channel_sets);
% 
%     % PCA map
%     subplot(nMaps + 1, length(phis), length(phis)*0 + nChannels_counter);
%     map = pca_map{nChannels_counter, reduced_dimensions}.score;
%     scatter(map(:, 1), plot_values, [], (1:size(channel_sets, 1)), '.');
%     title([num2str(nChannels) 'ch']);
%     ylabel(measure_string);
%     
%     % tSNE map
%     subplot(nMaps + 1, length(phis), length(phis)*1 + nChannels_counter);
%     map = tSNE_map{nChannels_counter, reduced_dimensions};
%     scatter(map(:, 1), plot_values, [], (1:size(channel_sets, 1)), '.');
%     ylabel(measure_string);
%     
%     % Path distance map
%     subplot(nMaps + 1, length(phis), length(phis)*2 + nChannels_counter);
%     centers = mean(channel_sets, 2); % Set centers - mean across channels in each set
%     distances = channel_set_distances(channel_sets); % Path distances
%     centers_offset_rand = (-0.1+(0.1--0.1).*rand(length(centers), 1)); % Values can vary by just under +- half of the minimum difference in centers
%     distances_offset_rand = (-0.4+(0.4--0.4).*rand(length(distances), 1)); % Values can vary by just under +- half of the minimum difference in distances
%     scatter(centers+centers_offset_rand, plot_values, [], (1:size(channel_sets, 1)), '.'); hold on;
%     
% %     % Channel set colorbar;
% %     subplot(nMaps + 1, length(phis), length(phis)*2 + nChannels_counter);
% %     scatter(pca_map{nChannels_counter, reduced_dimensions}.score(:, 1), tSNE_map{nChannels_counter, reduced_dimensions}(:, 1), [], (1:size(channel_sets, 1)), '.'); colorbar;
% %     title('Channel set mapping');
% %     xlabel('PC\_1'); ylabel('tSNE\_1D');
% %     axis('tight');
% end

%% Plot raw values against path 1D distance/centers

% flies = (1:13);
% 
% figure;
% %colormap('jet');
% condition_colours = 'rb';
% condition = 2;
% for nChannels_counter = 1 : length(phis)
%     nChannels = phis{nChannels_counter}.nChannels;
%     
%     if strcmp(measure_type, 'value')
%         
%         values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, flies, :, tau), 2), 3), [1 4 2 3]);
%         
%         plot_values = values(:, 1) - values(:, 2); % Awake - Anest
%         plot_values = values(:, 1)./values(:, 2);
%         plot_values = values(:, condition);
%         
%         value_sizes = (rescale(plot_values, 0.01, 200));
%         %value_sizes = rescale((plot_values).^5, 0.01, 400); % This changes the scale of values to be non-linear (in order to emphasise big/small values; assumes very small original values)
%         
%     else % strcmp(measure_type, 'class');
%         
%         values = mean(accuracies{nChannels_counter}.accuracies, 2);
%         values(values<58) = NaN;
%         
%         plot_values = values;
%         
%         value_sizes = rescale(plot_values, 0.01, 200);
%         
%     end
%     channel_sets = double(phis{nChannels_counter}.channel_sets);
%     
%     % Path distance map
%     subplot(2, length(phis), length(phis)*0 + nChannels_counter);
%     distances = channel_set_distances(channel_sets); % Path distances
%     distances_offset_rand = (-0.4+(0.4--0.4).*rand(length(distances), 1)); % Values can vary by just under +- half of the minimum difference in distances
%     scatter(distances+distances_offset_rand, plot_values, [], condition_colours(condition), '.'); hold on;
%     
%     % Path distance map
%     subplot(2, length(phis), length(phis)*1 + nChannels_counter);
%     centers = mean(channel_sets, 2); % Set centers - mean across channels in each set
%     centers_offset_rand = (-0.1+(0.1--0.1).*rand(length(centers), 1)); % Values can vary by just under +- half of the minimum difference in centers
%     scatter(centers+centers_offset_rand, plot_values, [], condition_colours(condition), '.'); hold on;
%     
% end

%% Plot change in values per channel
% Value for each channel = mean value across all networks which contain the channel

figure;
nChannels_colours = 'rgb';
for nChannels_counter = 1 : length(accuracies)
    
    if strcmp(measure_type, 'value')
        values = permute(mean(phis{nChannels_counter}.channel_means(:, :, :, :, tau), 2), [1 3 4 2]); % Average across trials
        
        plot_values = mean(values(:, :, 1) - values(:, :, 2), 2); title_string = ['Mean \Delta\Phi across sets containing channel X']; ylabel_string = [measure_string '_W - ' measure_string '_A'];
        axis_lims = [0 16 0 0.03]; % Awake - Anest
        
        plot_values = mean(values(:, :, 1)./values(:, :, 2), 2); plot_values_err = std(values(:, :, 2)./values(:, :, 1), [], 2) / sqrt(size(values, 2)); title_string = [measure_string '_W/' measure_string '_N across sets containing channel X']; ylabel_string = [measure_string '_W / ' measure_string '_N'];
        axis_lims = [0 16 1.6 4.8]; % [1 3.15] for phi3 within, [1.6 4.8] for phi3 across
        axis_lims = [0 16 1 2.5]; % for phi-star within
        axis_lims = [0 16 1 3.5]; % for phi-star across
        
    else % strcmp(measure_type, 'class')
        values = accuracies{nChannels_counter}.channel_means;
        plot_values = mean(values, 2); plot_values_err = std(values, [], 2) / sqrt(size(values, 2));
        title_string = ['Mean accuracy across sets containing channel X'];
        ylabel_string = '%';
        axis_lims = [0 16 50 75];
        axis_lims = [0 16 45 65];
    end
    
    % Plot
    errorbar((1:length(plot_values)), plot_values, plot_values_err./2, nChannels_colours(nChannels_counter)); hold on;
    %errorbar((1:length(plot_values)), plot_values, plot_values_err./2, [nChannels_colours(nChannels_counter) '.']); hold on;
    %errorbar((1:length(values_collapsed)), values_collapsed, values_collapsed_stderr);
    axis(axis_lims);
    legend('2ch', '3ch', '4ch', 'Location', 'northeastoutside');
    title(title_string);
    ylabel(ylabel_string); xlabel('channel');
end

%% Plot awake/anest comparison per channel (not classification)

figure;
nChannels_colours = 'rgb';
condition_lines = '-.';
axis_lims = [0 16 0 0.05]; % For phi3
%axis_lims = [0 16 0 0.006]; % For phi-star within
%axis_lims = [0 16 0 0.0007]; % For phi-star across
ylabel_string = measure_string;
title_string = [measure_string ' across sets containing channel X'];
for nChannels_counter = 1 : length(accuracies)
    
    values = permute(mean(phis{nChannels_counter}.channel_means(:, :, :, :, tau), 2), [1 3 4 2]); % Average across trials
    
    for condition = 1 : size(values, 3)
        plot_values = mean(values(:, :, condition), 2); plot_values_err = std(values(:, :, condition), [], 2) / sqrt(size(values, condition));
        errorbar((1:length(plot_values)), plot_values, plot_values_err./2, [nChannels_colours(nChannels_counter) condition_lines(condition)]); hold on;
    end
    
end

axis(axis_lims);
legend('2ch W', '3ch W', '4ch W', '2ch A', '3ch A', '4ch A', 'Location', 'northeastoutside');
title(title_string);
ylabel(ylabel_string); xlabel('channel');

%% Values for network vs means of individual values in network
% Plot value for each network against mean value of individual values in the network

figure;
colours = 'rb';
shapes = 'ox';
subplot_counter = 1;
for nChannels_counter = 2 : length(phis)
    subplot(1, 2, nChannels_counter - 1);
    
    nChannels = phis{nChannels_counter}.nChannels;
    channel_sets = phis{nChannels_counter}.channel_sets;
    
    for condition = 1 : 2
        if strcmp(measure_type, 'value')
            values = permute(mean(mean(phis{nChannels_counter}.phis(:, :, :, condition, tau), 2), 3), [1 4 2 3]);
            channel_means = permute(mean(mean(phis{1}.channel_means(:, :, :, condition), 2), 3), [1 4 2 3]);
        else % strcmp(measure_type, 'class')
            values = mean(accuracies{nChannels_counter}.accuracies, 2);
            channel_means = mean(accuracies{1}.channel_means, 2);
            channel_means = permute(mean(mean(phis{1}.channel_means(:, :, :, condition), 2), 3), [1 4 2 3]);
        end
        
        % Build values from average of single channel average values
        values_rebuilt = zeros(size(values));
        for channel_set = 1 : size(channel_sets, 1)
            channels = channel_sets(channel_set, :);
            values_rebuilt(channel_set) = mean(channel_means(channels));
        end
        
        % Plot
        scatter(values_rebuilt, values, [], colours(condition), shapes(condition)); hold on;
        
        % Stats
        [r, p] = corr(values_rebuilt, values);
        disp([num2str(nChannels_counter+1) 'Ch, condition ' num2str(condition) ': r=' num2str(r) ' p=' num2str(p)]);
        
    end
    
    title([num2str(nChannels) 'ch']);
    
    if strcmp(measure_type, 'value')
        legend('awake', 'anest', 'Location', 'northeastoutside');
        xlabel(['mean ' measure_string ' across channels']);
        ylabel(measure_string);
    else
        xlabel(['mean % across channels']);
        ylabel('%');
    end
    
end
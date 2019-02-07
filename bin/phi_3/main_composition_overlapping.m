%% DESCRIPTION

%{

Plots composition Hasse graph (phi-3 is z-axis)

Average across all flies, for non-global TPM

See figures/videos from http://www.eneuro.org/content/4/5/ENEURO.0085-17.2017

%}

%% Setup

output_file = 'animations/composition_overlapping_unpartitioned_setParamMap';
%'composition_overlapping_unpartitioned_orderMap';

flies = (1:13);

marker_size = 100;
nChannels = 4;

addpath('../');

%% Load

tic;
load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_phithree_nChannels4_globalTPM0.mat');
channel_sets = double(phis{1}.channel_sets);
toc

%% Get state-weighted compositions for all parameters

composition_phis = phis{1}.big_mips;

% Weight by state occurences (multiply phi by number of times the state occurred)
for partitioned = 1 : 2
    for concept = 1 : 15
        composition_phis(:, partitioned, concept, :, :, :, :) = ...
            permute(composition_phis(:, partitioned, concept, :, :, :, :), [1 4 5 6 7 2 3]) .* ...
            double(phis{1}.state_counters);
    end
end

% Sum across states
composition_phis = permute(sum(composition_phis, 1), [2 3 4 5 6 7 1]);

% Divide by total number of states (for weighted average)
% Assumes equal number of samples for all parameters
composition_phis = composition_phis ./ sum(phis{1}.state_counters(:, 1, 1, 1, 1));

% Unpartitioned - partitioned
composition_phis = permute(composition_phis(1, :, :, :, :, :) - composition_phis(2, :, :, :, :, :), [2 3 4 5 6 7 1]);
% Unpartitioned
%composition_phis = permute(composition_phis(1, :, :, :, :, :), [2 3 4 5 6 7 1]);
% Partitioned
%composition_phis = permute(composition_phis(2, :, :, :, :, :), [2 3 4 5 6 7 1]);
%% Setup Hasse graph

% Assumes x=1 is for highest order concept, x=2 is for second highest, etc

% Find x-y space limit requirements
nConcepts = zeros(size(channel_sets, 2), 1); % Number of points needed to show all concepts
channels = max(channel_sets(:));
for concept_order = 1 : size(channel_sets, 2)
    nConcepts(concept_order) = nchoosek(channels, concept_order);
end
concepts_cum = cumsum(nConcepts);
concepts_start = concepts_cum - nConcepts + 1;

% x-y space, with padding of 1
space_min = min(channel_sets(:))-1;
space_max = max(nConcepts)+1;

% x-y space here is based on concept-order and set-ID
x = zeros(1, sum(nConcepts));
y = zeros(size(x));
for concept_order = 1 : size(channel_sets, 2)
    x(concepts_start(concept_order):concepts_cum(concept_order)) = size(channel_sets, 2) - concept_order + 1; % highest order at x=1
    colours = x;
    spacing = linspace(space_min, space_max, nConcepts(concept_order)+2); % +2 is for padding of 1 on both sides
    y(concepts_start(concept_order):concepts_cum(concept_order)) = spacing(2:end-1);
end

% x-y space here is based on set centre and distance
x = zeros(1, sum(nConcepts));
y = zeros(size(x));
colours = zeros(size(x));
concept_counter = 1;
for concept_order = 1 : size(channel_sets, 2)
    concepts = nchoosek((1:channels), concept_order);
    for concept = 1 : length(concepts)
        x(concept_counter) = mean(concepts(concept, :)); % set center
        y(concept_counter) = channel_set_distance(concepts(concept, :)); % set distance
        concept_counter = concept_counter + 1;
    end
    colours(concepts_start(concept_order):concepts_cum(concept_order)) = size(channel_sets, 2) - concept_order + 1; % highest order at x=1
end

% Adjacency matrix
adj_mat = zeros(sum(nConcepts));
for concept_order = 1 : size(channel_sets, 2) - 1 % 4th order doesn't go anywhere
    from = nchoosek((1:channels), concept_order);
    to = nchoosek((1:channels), concept_order+1);
    for concept_lower = 1 : size(from, 1)
        for concept_higher = 1 : size(to, 1)
            if all(ismember(from(concept_lower, :), to(concept_higher, :)))
                adj_mat(...
                    concepts_start(concept_order) + concept_lower - 1,...
                    concepts_start(concept_order+1) + concept_higher - 1)...
                    = 1;
            end
        end
    end
end

figure;
scatter(x, y, [], colours, '.');

% g = digraph(adj_mat);
% plot(g, 'Layout', 'layered', 'ArrowSize', 0);

%% Convert phi composition into z axis

condition_titles = {'wake', 'anest'};

% Average across trials, flies (result is conditions x sets x concepts)
compositions = double(permute(mean(mean(composition_phis(:, :, :, flies, :), 3), 4), [5 2 1 3 4]));

% Find average concept phis across ALL sets which include the concept
comp_values = zeros(length(x), size(compositions, 1));
concept_counter = 1;
concept_displacement = 0; % For skipping lower order concepts
for concept_order = 1 : size(channel_sets, 2) % 4th-order concepts aren't shared
    concepts = nchoosek((min(channel_sets(:)):max(channel_sets(:))), concept_order);
    for concept = 1 : size(concepts, 1)
        value_sum = zeros(size(compositions, 1), 1); % values for each condition
        share_counter = 0;
        for network = 1 : size(channel_sets, 1)
            if all(ismember(concepts(concept, :), channel_sets(network, :))) % If concept is a subset of the channel set
                % Find phi of matching concept (concepts are ordered)
                network_concepts = nchoosek(channel_sets(network, :), concept_order);
                for network_concept = 1 : size(network_concepts, 1)
                    if all(ismember(concepts(concept, :), network_concepts(network_concept, :))) % If concept matches channel set concept
                        
%                         % Choose value based on max/min wake value
%                         if all(compositions(1, network, network_concept + concept_displacement) < value_sum(1))
%                             value_sum = compositions(:, network, network_concept + concept_displacement);
%                         end
                        
                        % Average values
                        value_sum = value_sum + compositions(:, network, network_concept + concept_displacement);
                        share_counter = share_counter + 1;
                        break;
                    end
                end
            end
        end
        
        % Selected value
%         comp_values(concept_counter, :) = value_sum;
        
        % Average values
        comp_values(concept_counter, :) = value_sum ./ share_counter;
        
        concept_counter = concept_counter + 1;
        
    end
    concept_displacement = concept_displacement + nchoosek(size(channel_sets, 2), concept_order);
end

%% Plot wake and anest on single axis

figure('pos', [0 0 1500 600]);
set(gcf, 'color', 'w');
condition_edges = [0 0 0; 0.9 0.9 0.9];
for condition = 1 : 2
    scatter3(x, y, comp_values(:, condition), marker_size, colours, 'o', 'filled', 'MarkerEdgeColor', condition_edges(condition, :));
    hold on;
    
%     % Draw lines
%     for concept = 1 : concepts_cum(end-1) % no lines from greatest-order concepts
%         destinations = find(adj_mat(concept, :));
%         for dest_counter = 1 : length(destinations)
%             dest = destinations(dest_counter);
%             tmp = line(...
%                 [x(concept) x(dest)],...
%                 [y(concept), y(dest)],...
%                 [comp_values(concept, condition), comp_values(dest, condition)],...
%                 'LineStyle', '-',...
%                 'Color', 'k',...
%                 'LineWidth', 0.5 ...
%                 );
%             tmp.Color(4) = 0.05;
%         end
%     end
    
end

zlabel('\phi');
xlabel('set-centre');
ylabel('set-distance');

set(gca, 'YTick', [min(y) max(y)], 'XTick', [min(x) max(x)], 'ZTick', linspace(min(comp_values(:)), max(comp_values(:)), 3));

%axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 min(comp_values(:)) max(comp_values(:))]);
axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 0 max(comp_values(:))]);

box off
grid off
%axis square
axis vis3d

title('Mean concept \phi across all channel sets');

%% Plot wake - anest


figure('pos', [0 0 1500 600]);
set(gcf, 'color', 'w');

scatter3(x, y, comp_values(:, 1) - comp_values(:, 2), marker_size, colours, '.');
hold on;

% Draw lines
for concept = 1 : concepts_cum(end-1) % no lines from greatest-order concepts
    destinations = find(adj_mat(concept, :));
    for dest_counter = 1 : length(destinations)
        dest = destinations(dest_counter);
        tmp = line(...
            [x(concept) x(dest)],...
            [y(concept), y(dest)],...
            [comp_values(concept, 1) - comp_values(concept, 2), comp_values(dest, 1) - comp_values(dest, 2)],...
            'LineStyle', '-',...
            'Color', 'k',...
            'LineWidth', 0.5 ...
            );
        tmp.Color(4) = 0.05;
    end
end

zlabel('\phi');
xlabel('x');
ylabel('y');
axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 min(comp_values(:,1)-comp_values(:,2)) max(comp_values(:,1)-comp_values(:,2))]);

%title([condition_titles{condition} ': \Phi=' num2str(phis{condition}.phi)]);
title('wake-anest');

set(gca, 'YTick', [min(y) max(y)], 'XTick', [min(x) max(x)], 'ZTick', linspace(0, max(compositions(:)), 3));
set(gca, 'YTickLabel', [], 'XTickLabel', []);

box off
grid off
%axis square
axis vis3d


%% Plot wake and anest

figure('pos', [0 0 1500 600]);
set(gcf, 'color', 'w');
subplots = zeros(1, 4);
subplot_counter = 1;
for condition = 1 : 2
    subplots(condition) = subplot(1, 2, subplot_counter);
    scatter3(x, y, comp_values(:, condition), marker_size, colours, '.');
    hold on;
    
    % Draw lines
    for concept = 1 : concepts_cum(end-1) % no lines from greatest-order concepts
        destinations = find(adj_mat(concept, :));
        for dest_counter = 1 : length(destinations)
            dest = destinations(dest_counter);
            tmp = line(...
                [x(concept) x(dest)],...
                [y(concept), y(dest)],...
                [comp_values(concept, condition), comp_values(dest, condition)],...
                'LineStyle', '-',...
                'Color', 'k',...
                'LineWidth', 0.5 ...
                );
            tmp.Color(4) = 0.05;
        end
    end
    
    zlabel('\phi');
    xlabel('set-centre');
    ylabel('set-distance');
    axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 min(comp_values(:)) max(comp_values(:))]);
    
    %title([condition_titles{condition} ': \Phi=' num2str(phis{condition}.phi)]);
    title([condition_titles{condition}]);
    
    set(gca, 'YTick', [min(y) max(y)], 'XTick', [min(x) max(x)], 'ZTick', linspace(min(comp_values(:)), max(comp_values(:)), 3));
    %set(gca, 'YTickLabel', [], 'XTickLabel', []);
    
    box off
    grid off
    %axis square
    axis vis3d
    subplot_counter = subplot_counter + 1;
end

linkprop(subplots, {'CameraPosition','CameraUpVector'});

%% Rotate and turn into video frames

% Angles of elevation for top and side views
el_top = 89; % at 90-89, the y axis location flips from one side of the plot to the other, so start from 89
el_hor = 0;
az_min = 0;
az_max = 360 + az_min;


% Rotation angles
top2hor = linspace(el_top, el_hor, 100); % rotate from top view to side view
hor2hor = linspace(az_min, az_max, 360); % rotate side view
hor2top = linspace(el_hor, el_top, 100); % rotate from side view to top view

top2hor = linspace(el_top, el_hor, 100/4); % rotate from top view to side view
hor2hor = linspace(az_min, az_max, 360/4); % rotate side view
hor2top = linspace(el_hor, el_top, 100/4); % rotate from side view to top view

% Pad vertical and horizontal rotations
azimuth = [zeros(size(top2hor))+az_min hor2hor zeros(size(hor2top))+az_min]; % rotation during side view
elevation = [top2hor zeros(size(hor2hor))+el_hor hor2top]; % rotation from top to side and back

frames = cell(size(azimuth));
for view_angle = 1 : length(azimuth)
    disp(view_angle);
    tic;
    view([azimuth(view_angle) elevation(view_angle)]);
    drawnow
    frame = getframe(gcf);
    frames{view_angle} = frame2im(frame);
    toc
end

%% Write frames into gif

video_duration = 20; % in seconds
frame_duration = video_duration / length(azimuth);

for frame = 1 : length(frames)
    [mapped_frame, map] = rgb2ind(frames{frame}, 256);
    if frame == 1
        imwrite(mapped_frame, map, [output_file '.gif'], 'gif', 'LoopCount', Inf, 'DelayTime', frame_duration);
    else
        imwrite(mapped_frame, map, [output_file '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', frame_duration);
    end
end

%% Write frames into video

% create the video writer with 1 fps
writerObj = VideoWriter([output_file '.mp4'], 'MPEG-4');
writerObj.FrameRate = 1 / (video_duration / length(frames));

% open the video writer
open(writerObj);

% write the frames to the video
for u=1:length(frames)
    % convert the image to a frame
    frame = im2frame(frames{u});
    
    % write to video
    writeVideo(writerObj, frame);
end

% close the writer object
close(writerObj);
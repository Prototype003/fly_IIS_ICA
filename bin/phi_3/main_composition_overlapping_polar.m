%% DESCRIPTION

%{

Plots composition Hasse graph (phi-3 is z-axis)

Average across all flies, for non-global TPM

See figures/videos from http://www.eneuro.org/content/4/5/ENEURO.0085-17.2017

%}

%% Setup

output_file = 'animations/composition_overlapping_unpartitioned_setParamMap_polar';
%'composition_overlapping_unpartitioned_orderMap';

flies = (1:13);
trials = (1:8);

marker_size = 10;
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
    xspacing = linspace(space_min, space_max, nConcepts(concept_order)+2); % +2 is for padding of 1 on both sides
    y(concepts_start(concept_order):concepts_cum(concept_order)) = xspacing(2:end-1);
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

% figure;
% scatter(x, y, [], colours, '.');

% x = rescale(x, 0, 360);
% xp = (y+1) .* cosd(x);
% yp = (y+1) .* sind(x);

y = rescale(y, 0, 360);
xp = (max(x)-x+1) .* cosd(y); % max(x) - x reverses set center, so most central becomes most peripheral
yp = (max(x)-x+1) .* sind(y);

x = xp; y = yp;
colours = 4-colours; % Invert the colour order
% 
% figure;
% scatter(xp, yp, [], colours, '.');

% g = digraph(adj_mat);
% plot(g, 'Layout', 'layered', 'ArrowSize', 0);

%% Convert phi composition into z axis

condition_titles = {'wake', 'anest'};

% Average across trials, flies (result is conditions x sets x concepts)
compositions = double(permute(mean(mean(composition_phis(:, :, trials, flies, :), 3), 4), [5 2 1 3 4]));

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
%                         if all(compositions(1, network, network_concept + concept_displacement) > value_sum(1))
%                             value_sum = compositions(:, network, network_concept + concept_displacement);
%                         end
                        
                        % Average values
                        value_sum = value_sum + compositions(:, network, network_concept + concept_displacement);
                        share_counter = share_counter + 1;
                        
                        break; % The concept only occurs once
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

figure('pos', [0 0 800 600]);
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
    
    zlabel('\phi');
    xlabel('set-centre');
    ylabel('set-distance');
    
    %title([condition_titles{condition} ': \Phi=' num2str(phis{condition}.phi)]);
    %title([condition_titles{condition}]);
    title('Mean concept \phi across all channel sets');
    
    %set(gca, 'YTick', [min(y) max(y)], 'XTick', [min(x) max(x)], 'ZTick', linspace(min(comp_values(:)), max(comp_values(:)), 3));
    %set(gca, 'YTickLabel', [], 'XTickLabel', []);
    
end

%set(gca, 'YTick', [min(y) max(y)], 'XTick', [min(x) max(x)], 'ZTick', linspace(min(comp_values(:)), max(comp_values(:)), 3));
set(gca, 'ZTick', linspace(0, max(comp_values(:)), 3));

%axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 min(comp_values(:)) max(comp_values(:))]);
axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 0 max(comp_values(:))]);

box off
grid off
%axis square
axis vis3d

title('Mean concept \phi across all channel sets');

handle = gca;
handle.XAxis.FirstCrossoverValue = 0; % Doesn't work for low elevation
handle.YAxis.FirstCrossoverValue = 0; % Doesn't work for low elevation
handle.ZAxis.FirstCrossoverValue = min(x);
handle.ZAxis.SecondCrossoverValue = min(y);
handle.ZAxis.TickLength = [0 0];
handle.ZAxis.Visible = 'off';

% Centre the plot at (x=y=0)
xlimits = xlim; xlim([-max(xlimits) max(xlimits)]);
ylimits = ylim; ylim([-max(ylimits) max(ylimits)]);
ztick_length = 0.5; % Should make dynamic based on range of x/y axis
xspacing = 3; % Should make dynamic based on range of x/y axis

% Manually draw x/y axes (crossover value doesn't work at low elevation)
set(gca, 'XColor', 'none', 'YColor', 'none');
line([0 0], [0 0], [0 1], 'Color', 'k'); % z
for zval = linspace(0, max(comp_values(:)), 3)
    line([0 -0.5], [0 0], [zval zval], 'Color', 'k');
end
line(xlimits, [0 0], [0 0], 'Color', 'k');
line([0 0], ylimits, [0 0], 'Color', 'k');

% Custom x-axis ticks
xtick_length = 0.00015; % Should be made dynamic based on range of z
zspacing = 0.00025;
zlimits = zlim; zlim([zlimits(1)-xtick_length zlimits(2)]); % Extend by length of tick mark
zlimits = zlim;
xvals = [1 15];
xlabels = {'peripheral', 'central'};
for counter = 1 : length(xvals)
    xval = xvals(counter);
    line([xval xval], [0 0], [0 zlimits(1)], 'Color', 'k');
    label = text(max(xvals)-xval+1, 0, zlimits(1)-zspacing, xlabels{counter}, 'HorizontalAlignment', 'center');
end

% Custom z-axis
xlimits = xlim; % Wider than previously (using new limits)
ylimits = ylim; % Wider than previously (using new limits)
line([xlimits(1) xlimits(1)]+ztick_length, [0 0], [0 1], 'Color', 'k');
for zval = linspace(0, max(comp_values(:)), 3)
    line([xlimits(1) xlimits(1)-ztick_length]+ztick_length, [0 0], [zval zval], 'Color', 'k');
end
tick_points = linspace(0, max(comp_values(:)), 3);
label = text(xlimits(1)-xspacing*3, 0, tick_points(2), '\phi', 'HorizontalAlignment', 'center');
label.Rotation = 90;
for zval = tick_points
    label = text(xlimits(1)-xspacing, 0, zval, num2str(zval, '%0.4f'), 'HorizontalAlignment', 'center');
end

legend('wake', 'anest', 'Location', 'southeast');

%% Plot wake - anest

figure('pos', [0 0 1500 600]);
set(gcf, 'color', 'w');

scatter3(x, y, comp_values(:, 1) - comp_values(:, 2), marker_size, colours, 'o', 'filled', 'MarkerEdgeColor', [0 0 0]);
hold on;

% % Draw lines
% for concept = 1 : concepts_cum(end-1) % no lines from greatest-order concepts
%     destinations = find(adj_mat(concept, :));
%     for dest_counter = 1 : length(destinations)
%         dest = destinations(dest_counter);
%         tmp = line(...
%             [x(concept) x(dest)],...
%             [y(concept), y(dest)],...
%             [comp_values(concept, 1) - comp_values(concept, 2), comp_values(dest, 1) - comp_values(dest, 2)],...
%             'LineStyle', '-',...
%             'Color', 'k',...
%             'LineWidth', 0.5 ...
%             );
%         tmp.Color(4) = 0.05;
%     end
% end

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

%% Rotate and turn into video frames

% Angles of elevation for top and side views
el_top = 90; % at 90-89, the y axis location flips from one side of the plot to the other, so start from 89
el_hor = 5;
az_min = 0;
az_max = 360 + az_min;


% Rotation angles
top2hor = linspace(el_top, el_hor, 100); % rotate from top view to side view
hor2hor = linspace(az_min, az_max, 360); % rotate side view
hor2top = linspace(el_hor, el_top, 100); % rotate from side view to top view

top2hor = linspace(el_top, el_hor, 100); % rotate from top view to side view
hor2hor = linspace(az_min, az_max, 360); % rotate side view
hor2top = linspace(el_hor, el_top, 100); % rotate from side view to top view

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
%% DESCRIPTION

%{

Plots composition Hasse graph (phi-3 is z-axis)

See figures/videos from http://www.eneuro.org/content/4/5/ENEURO.0085-17.2017

%}

%% Setup

output_file = '../figures/iis_example';

marker_size = 500;

%% Load phi compositions

load('results/split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_ICA_phithree_nChannels4_globalTPM0.mat');

%% Average across states, trials

composition_phis = phis{1}.big_mips;
%nonzero = composition_phis(composition_phis~=0);
%composition_phis = log(composition_phis+min(nonzero));
%composition_phis = log(composition_phis+1);

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

% Average across trials
composition_phis = mean(composition_phis, 4);

% concepts x channel-sets x trials x flies x conditions x part-types
dims = size(composition_phis);
compositions = zeros([dims(2:end) 3]);

% Unpartitioned
composition_unpart = permute(composition_phis(1, :, :, :, :, :), [2 3 4 5 6 7 1]);
compositions(:, :, :, :, :, 1) = composition_unpart;
% Partitioned
composition_part = permute(composition_phis(2, :, :, :, :, :), [2 3 4 5 6 7 1]);
compositions(:, :, :, :, :, 2) = composition_part;
% Unpartitioned - partitioned
composition_diff = permute(composition_phis(1, :, :, :, :, :) - composition_phis(2, :, :, :, :, :), [2 3 4 5 6 7 1]);
compositions(:, :, :, :, :, 3) = composition_diff;

%% Setup Hasse graph

% Hard coded coordinates assumes order:
% ABCD BCD ACD ABD ABC CD BD BC AD AC AB D C B A

% x-y space, with padding of 1
space_min = 0;
space_max = 7;

% Each row corresponds to a mechanism size
order_colours = [233 163 50; 131 197 90; 70 146 207; 69 60 151] / 255;

colours = [1 2 2 2 2 3 3 3 3 3 3 4 4 4 4];
colours = cat(1, order_colours(4, :),...
    repmat(order_colours(3, :), [4 1]),...
    repmat(order_colours(2, :), [6 1]),...
    repmat(order_colours(1, :), [4 1])); 

space_1 = space_max / 2;
space_4 = linspace(space_min, space_max, 4+2); space_4 = space_4(2:end-1);
space_6 = linspace(space_min, space_max, 6+2); space_6 = space_6(2:end-1);

y = [space_1 space_4 space_6 space_4];
x = [1 zeros(size(space_4))+2 zeros(size(space_6))+3 zeros(size(space_4))+4];

% figure;
% scatter(x, y);

%% Convert phi composition into z axis

nChannels = 4;

% Select channel set, fly to plot
channel_set = 1;
fly = 11;
part_type = 1;
compositions_plot = permute(compositions(:, channel_set, :, fly, :, part_type), [5 1 2 3 4 6]);

% Repeat the first composition
compositions_plot = cat(1, compositions_plot(1, :), compositions_plot);

% Set the first composition to all min values
vals = compositions_plot(2:end, :);
compositions_plot(1, :) = min(vals(:));

% Higher orders first (because lines are in reverse order)
compositions_plot = fliplr(compositions_plot);

% Hardcoded lines for Hasse Diagram
% Cleanest way to plot these?
% Cell array for each concept, contains list of indexes to which a line should be drawn
% ABCD to BCD, ACD, ABD, ABC (1 > 2, 3, 4, 5)
% BCD to CD, BD, BC (2 > 6, 7, 8)
% ACD to CD, AD, AC (3 > 6, 9, 10)
% ABD to BD, AD, AB (4 > 7, 9, 11)
% ABC to BC, AC, AB (5 > 8, 10, 11)
% CD to D, C (6 > 12, 13)
% BD to D, B (7 > 12, 14)
% BC to C, B (8 > 13, 14)
% AD to D, A (9 > 12, 15)
% AC to C, A (10 > 13, 15)
% AB to B, A (11 > 14, 15)

lines = cell(length(y) - nChannels, 1);
lines{1} = [2 3 4 5];
lines{2} = [6 7 8];
lines{3} = [6 9 10];
lines{4} = [7 9 11];
lines{5} = [8 10 11];
lines{6} = [12 13];
lines{7} = [12 14];
lines{8} = [13 14];
lines{9} = [12 15];
lines{10} = [13 15];
lines{11} = [14 15];

figure;
subplots = zeros(1, size(compositions_plot, 1));
subplot_counter = 1;
subtitles = {'', 'wake', 'anest'};
for condition = 1 : size(compositions_plot, 1)
    subplots(condition) = subplot(1, size(compositions_plot, 1), subplot_counter);
    scatter3(x, y, compositions_plot(condition, :), marker_size, colours, '.');
    
    % Draw lines
    for source = 1:length(lines)
        for dest = lines{source}
            line([x(source) x(dest)], [y(source) y(dest)], [compositions_plot(condition, source) compositions_plot(condition, dest)], 'Color', 'k');
        end
    end
    
    zlabel('\phi', 'Rotation', 90);
    ylabel('y');
    if subplot_counter == 1
        xlabel('mechanism size');
    end
    
    % Set axis limits (based on other plots)
    axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 min(compositions_plot(:)) max(compositions_plot(:))]);
    set(gca, 'ZTick', linspace(0, max(compositions_plot(:)), 3));
    
%     % Set axis limits (based on own values)
%     axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 min(compositions(:)) max(compositions(condition, :))]);
%     set(gca, 'ZTick', linspace(0, max(compositions(end, :)), 3));
    
%     if condition <=2
%         title([condition_titles{condition} ': \Phi=' num2str(phis{condition}.phi)]);
%     end
    
    title(subtitles{condition});
    
    set(gca, 'YTick', [min(y) max(y)], 'XTick', [min(x) max(x)]);
    set(gca, 'YTickLabel', [], 'XTickLabel', [4 1]);
    
    box off
    grid on
    axis square
    axis vis3d
    view([10+135 10]);
    if condition == 1
        view([180 90]);
    end
    subplot_counter = subplot_counter + 1;
end

linkprop(subplots(2:end), {'CameraPosition','CameraUpVector'});

%% Output still figure

figure_name = '../figures/fig2_raw';

set(gcf, 'PaperOrientation', 'Landscape');

print(figure_name, '-dsvg', '-painters'); % SVG
print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
print(figure_name, '-dpng'); % PNG

%% Rotate and turn into video frames

% Angles of elevation for top and side views
el_top = 89; % at 90-89, the y axis location flips from one side of the plot to the other, so start from 89
el_hor = 0;
az_min = 180;
az_max = 360 + az_min;

% Rotation angles
top2hor = linspace(el_top, el_hor, 100); % rotate from top view to side view
hor2hor = linspace(az_min, az_max, 360); % rotate side view
hor2top = linspace(el_hor, el_top, 100); % rotate from side view to top view

% Pad vertical and horizontal rotations
azimuth = [zeros(size(top2hor))+az_min hor2hor zeros(size(hor2top))+az_min]; % rotation during side view
elevation = [top2hor zeros(size(hor2hor))+el_hor hor2top]; % rotation from top to side and back

frames = cell(size(azimuth));
for view_angle = 1 : length(azimuth)
    view([azimuth(view_angle) elevation(view_angle)]);
    drawnow
    frame = getframe(gcf);
    frames{view_angle} = frame2im(frame);
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
 writerObj = VideoWriter([output_file], 'MPEG-4');
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
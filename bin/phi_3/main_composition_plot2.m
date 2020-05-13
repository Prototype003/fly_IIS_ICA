%% DESCRIPTION

%{

Plots composition Hasse graph (phi-3 is z-axis)

Difference from original - draws lines to reflect MIPs, instead of all
related mechanisms from the previous order

See figures/videos from http://www.eneuro.org/content/4/5/ENEURO.0085-17.2017

%}

%% Setup

output_file = 'composition_both';

marker_size = 500;

%% Setup Hasse graph

% Hard coded coordinates assumes order:
% ABCD BCD ACD ABD ABC CD BD BC AD AC AB D C B A

% x-y space, with padding of 1
space_min = 0;
space_max = 7;

colours = [1 2 2 2 2 3 3 3 3 3 3 4 4 4 4];

space_1 = space_max / 2;
space_4 = linspace(space_min, space_max, 4+2); space_4 = space_4(2:end-1);
space_6 = linspace(space_min, space_max, 6+2); space_6 = space_6(2:end-1);

y = [space_1 space_4 space_6 space_4];
x = [1 zeros(size(space_4))+2 zeros(size(space_6))+3 zeros(size(space_4))+4];

% figure;
% scatter(x, y);

%% Convert phi composition into z axis

load('split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f01c1tau4s0876t1_example.mat');


compositions = zeros(4, length(mech_list)-1);

% Get mechanism phi values
%for m = 2 : length(mech_list) % first mechanism is a null mechanism
%end

%% Convert phi composition into z axis

condition_titles = {'wake', 'anest'};
phis = cell(1, 2);

nChannels = 4;

% Nearest mean big-phi? - set 129, 84.3195%
% SVM big-phi? - set 876, 81.95%
% SVM 4th order concept? - set 74, 87.87%
% SVM full composition? - set 234, 84.32%
channel_set = 876;

file_prefix = 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels';

[phis{1}, wake_d, wake_u, wake_p] = composition_table(nChannels, [file_prefix num2str(nChannels) '_globalTPM0_f01c1tau4s' sprintf('%04d', channel_set) 't1.mat']);
[phis{2}, anest_d, anest_u, anest_p] = composition_table(nChannels, [file_prefix num2str(nChannels) '_globalTPM0_f01c1tau4s' sprintf('%04d', channel_set) 't1.mat']);

%[phis{1}, wake] = composition_table(nChannels, 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f01c1tau4tauOffset0s0002t1.mat');
%[phis{2}, anest] = composition_table(nChannels, 'split2250_bipolarRerefType1_lineNoiseRemoved_postPuffpreStim_nChannels4_globalTPM0_f01c2tau4tauOffset0s0002t1.mat');

%compositions = [wake_u(end, :); anest_u(end, :); wake_p(end, :); anest_p(end, :)];
compositions = [zeros(1, size(wake_u, 2)); wake_u(end, :); wake_p(end, :); wake_d(end, :)];
compositions = fliplr(compositions);

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
subplots = zeros(1, size(compositions, 1));
subplot_counter = 1;
subtitles = {'', 'unpart.', 'part.', 'unpart. - part.'};
for condition = 1 : size(compositions, 1)
    subplots(condition) = subplot(1, size(compositions, 1), subplot_counter);
    scatter3(x, y, compositions(condition, :), marker_size, colours, '.');
    
    % Draw lines
    for source = 1:length(lines)
        for dest = lines{source}
            line([x(source) x(dest)], [y(source) y(dest)], [compositions(condition, source) compositions(condition, dest)], 'Color', 'k');
        end
    end
    
    zlabel('\phi', 'Rotation', 90);
    xlabel('x');
    ylabel('y');
    
    if subplot_counter < size(compositions, 1)
        axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 0 max(compositions(:))]);
        set(gca, 'ZTick', linspace(0, max(compositions(:)), 3));
    else
        axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1 0 max(compositions(condition, :))]);
        set(gca, 'ZTick', linspace(0, max(compositions(end, :)), 3));
    end
%     if condition <=2
%         title([condition_titles{condition} ': \Phi=' num2str(phis{condition}.phi)]);
%     end
    
    title(subtitles{condition});
    
    set(gca, 'YTick', [min(y) max(y)], 'XTick', [min(x) max(x)]);
    set(gca, 'YTickLabel', [], 'XTickLabel', []);
    
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

%linkprop(subplots, {'CameraPosition','CameraUpVector'});

%% Output still figure

% figure_name = 'figures/fig2_raw';
% 
% set(gcf, 'PaperOrientation', 'Landscape');
% 
% print(figure_name, '-dsvg', '-painters'); % SVG
% print(figure_name, '-dpdf', '-painters', '-bestfit'); % PDF
% print(figure_name, '-dpng'); % PNG

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
 writerObj = VideoWriter([output_file '.avi']);
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
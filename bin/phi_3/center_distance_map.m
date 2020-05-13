function [values_map_interp, values_map, centers, distances, centers_axis, distances_axis] = center_distance_map(channel_sets, values)
% Maps sets of linear array channels according to 2 dimensions - center of
% the channels and sum of pairwise distances between all pairs of channels
% in the set
%
% Inputs:
%   channel_sets = matrix (number of channel sets x number of channels in
%       each set)
%   values = vector (number of channel set x 1) giving values associated
%       with each channel set
%
% Outputs:
%   values_map_interp = matrix (distances x centers);
%       granularity/resolution of distances and centers is determined by
%       number of channels in each set; 'gaps' in the mapping are
%       interpolated linearly
%   values_map = matrix (distances x centers); values in the 2D mapping,
%       before interpolation; overlapping values at the same mapping are
%       averaged (0 when there are no values in the mapped coordinate)
%   centers = vector (number of channel sets x 1) giving the centers of
%       each channel set
%   distances = vector (number of channels sets x 1) giving the total
%       pairwise distances of each channel set

% Set center / set path distance mapping
centers = mean(channel_sets, 2); % Mean across channels in each set
distances = channel_set_distances(channel_sets);

% Store mapping
values_mapped.centers = centers;
values_mapped.distances = distances;
values_mapped.values = values;

% Find minimum delta among map values
centers_deltas = pdist2(unique(centers), unique(centers));
centers_deltas(centers_deltas==0) = NaN; % remove 0 distances
centers_delta = min(centers_deltas(:));
distances_deltas = pdist2(unique(distances), unique(distances));
distances_deltas(distances_deltas==0) = NaN; % remove 0 distances
distances_delta = min(distances_deltas(:));

% Create space mapping (used for indexing into the mapped space)
%     centers_axis = (min(centers) - (2*centers_delta) : centers_delta : max(centers) + (2*centers_delta)); % - and + delta is for padding
%     distances_axis = (min(distances) - (2*distances_delta) : distances_delta : max(distances) + (2*distances_delta)); % - and + delta is for padding
%     [centers_map, distances_map] = meshgrid(centers_axis, distances_axis);
centers_axis = (min(channel_sets(:)) : centers_delta : max(channel_sets(:))+centers_delta); % range is from first channel to last channel; +delta is for padding (because of how pcolor works)
distances_axis = (min(distances) : distances_delta : max(distances)+distances_delta); % range is from min to max
[centers_map, distances_map] = meshgrid(centers_axis, distances_axis);

% Map channel sets to 2D matrix (and average overlaps)
% Create mapped space
values_map = zeros(size(centers_map)); % Will sum all values with the same coordinates
values_map_counter = zeros(size(centers_map)); % Keeps count in each coordinate as to how many values have that coordinate

% Populate mapped space
for value_counter = 1 : length(values)
    x = find(abs(centers_axis - centers(value_counter)) < 0.00001, 1); % This gives the mapped x location
    y = find(distances_axis == distances(value_counter), 1); % This gives the mapped y location
    values_map(y, x) = values_map(y, x) + values(value_counter); % matrix is (rows, columns), corresponding to (y, x)
    values_map_counter(y, x) = values_map_counter(y, x) + 1;
end

values_map_interp = values_map ./ values_map_counter; % We will use the NaN values for the black background (and maybe interpolation)

% Interpolate values (linearly - each NaN will turn into the average of the cells immediately adjacent to it, excluding diagonals)
min_value = min(values_map_interp(:)); % We will remove interpolations which are less than the original min value in the plot (to avoid 'ghost/shadow interpolations');
values_adjacent = zeros(4, 1);
for y = 2 : size(values_map_interp, 1) - 1
    for x = 2 : size(values_map_interp, 2) - 1
        if isnan(values_map_interp(y, x))
            values_adjacent(1) = values_map_interp(y, x-1);
            values_adjacent(2) = values_map_interp(y, x+1);
            values_adjacent(3) = values_map_interp(y-1, x);
            values_adjacent(4) = values_map_interp(y+1, x);
            if sum(~isnan(values_adjacent)) > 2
                %values_adjacent(isnan(values_adjacent)) = 0; % Deal with NaNs
                %values_adjacent(values_adjacent < min_value) = 0; % Don't include 'shadow' interpolations
                interpolated = mean(values_adjacent(~isnan(values_adjacent)));
                if interpolated == 0
                    values_map_interp(y, x) = NaN; % Place NaN instead of 0 to avoid affecting the color scaling
                else
                    values_map_interp(y, x) = interpolated;
                end
            end
        end
    end
end

values_map = values_map ./ values_map_counter; % average
values_map(isnan(values_map)) = 0; % replace NaNs with 0

centers_axis = centers_map(1, :);
distances_axis = distances_map(:, 1);

end


function fig = plot_map(fig, x, path)
%PLOT_MAP
%
% fig = plot_map([], map, path)
% fig = plot_map(fig, map, path)


if nargin < 3 || isempty(path)
    path = [];
end
if isa(x, 'VoxelMap')
    map = x;
else
    map = [];
end

%% Get figure/axes limits.
if ~isempty(path)
    margin = [40 40 3]';
    x_lims = [min(path, [], 2)-(margin+[0 0 1]'), max(path, [], 2) + margin];
elseif ~isempty(map)
    % Use occupied voxels to set limits.
    [x, ~] = map.get_voxels();
    x_lims = [min(x, [], 2) max(x, [], 2)];
    
else
    x_lims = [min(x, [], 2) max(x, [], 2)];
end

%% Initialize figure.
if isempty(fig)
    fig = figure('Name', 'Map');
else
    figure(fig);
end
hold on;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;
axis off;
xlim(x_lims(1,:));
ylim(x_lims(2,:));
zlim(x_lims(3,:));


camlight('right')
lighting gouraud;

%% Plot map points/surface.

colors = parula(256);

if ~isempty(map)
%     threshold = map.occupied_threshold;
    threshold = 0;
    colormap(colors);
    if isempty(path)
        % Use all points within limits (may be slow).
        [gx, gy, gz] = meshgrid(...
            x_lims(1,1):map.voxel_size:x_lims(1,2), ...
            x_lims(2,1):map.voxel_size:x_lims(2,2), ...
            x_lims(3,1):map.voxel_size:x_lims(3,2));
        x = [gx(:)'; gy(:)'; gz(:)'];
        [~, val] = map.get_voxels(x);
%         val(isnan(val)) = min(val);
        val(isnan(val)) = 0;
        val = reshape(val, size(gx));
        norm_z = clip(map_intervals(x(3,:)', x_lims(3,:), [0 1]), [0 1]);
        figure(fig);
        isosurface(gx, gy, gz, val, threshold, reshape(norm_z, size(gx)), 'noshare');
    else
        % Use neighborhood along path, skip position too close to each other.
        for i = 1:size(path, 2)
            if i > 1 && i < size(path, 2) && all(abs(path(:,i+1) - path(:,last_i)) < margin/2)
                continue;
            end
            last_i = i;
            [gx, gy, gz] = meshgrid(...
            path(1,i)-margin(1):map.voxel_size:path(1,i)+margin(1), ...
            path(2,i)-margin(2):map.voxel_size:path(2,i)+margin(2), ...
            path(3,i)-margin(3):map.voxel_size:path(3,i)+margin(3));
            x = [gx(:)'; gy(:)'; gz(:)'];
            [~, val] = map.get_voxels(x);
            val(isnan(val)) = 0;
            val = reshape(val, size(gx));
            norm_z = clip(map_intervals(x(3,:)', x_lims(3,:), [0 1]), [0 1]);
            figure(fig);
            isosurface(gx, gy, gz, val, threshold, reshape(norm_z, size(gx)),'noshare');
        end
    end
else
    norm_z = clip(map_intervals(x(3,:)', x_lims(3,:), [0 1]), [0 1]);
    figure(fig);
    plot3_color(x', map_rgb(norm_z, colors));
end

if ~isempty(path)
     plot3(path(1, :), path(2, :), path(3, :), 'ko-', 'LineWidth', 1);
     plot3(path(1, end), path(2, end), path(3, end), 'ro-', 'LineWidth', 1);

end

    function y = map_intervals(x, from, to)
        y = (x - from(1)) / diff(from);
        y = diff(to) * y + to(1);
    end

    function y = clip(x, lims)
        y = max(min(x, lims(2)), lims(1));
    end

    function rgb_z = map_rgb(z, rgb)
        z = floor(max(min(z, 1-eps), 0) * size(rgb, 1)) + 1;
        rgb_z = rgb(z, :);
    end



end

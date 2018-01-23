function map = build_voxel_map(map, base_dir, T, range)
%BUILD_VOXEL_MAP
%
% map = build_voxel_map(map, base_dir, T)
% map = build_voxel_map(map, base_dir, T, range)
%
% Input:
% - map       An initial VoxelMap to be updated, or [].
% - base_dir  Base data dir containing
%             [base_dir '/velodyne_points/data/%010i.bin'].
% - T         Cell arrya of 4-by-4 transformation matrices.
% - range     Minimum and maximum measurement range to be used.
%
% Output:
% - map       The updated VoxelMap instance.
%

if nargin < 4 || isempty(range)
    range = [0 inf];
end
assert(isa(range, 'double'));
assert(numel(range) == 2);
if isempty(map)
    map = VoxelMap();
end

for i = 1:numel(T)
    if isempty(T{i})
        continue;
    end
    cloud_name = sprintf('%010i.bin', i-1);
    cloud_path = [base_dir '/velodyne_points/data/' cloud_name];
    x = read_velo(cloud_path, range);
    % Update map with transformed points.
    x_map = p2e(T{i} * e2p(x));
    x_orig = p2e(T{i} * [0 0 0 1]');
    t = tic();
    map.update_lines(x_orig, x_map);
    fprintf('Map update with %i points from %s: %.3f s.\n', size(x_map, 2), cloud_name, toc(t));
end

end

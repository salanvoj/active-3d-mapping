function [map_merged, map_conf, map_meas, map_gt, map_measurable, path] = run_pipeline(map_dir, net_path, planning, draw)

% Voxel map params
voxel_size = 0.2;
free_update = -1;
hit_update = 2;
occupied_threshold = 2;
measurement_range = [2.5 100];

% Load transforms from velo to map, create sensor path.
velo_to_map.T = load_velo_to_map(map_dir);
path = cell2mat(cellfun(@(T) T(1:3,4), velo_to_map.T, 'UniformOutput', false));

% Build ground-truth map.
if exist(strcat(map_dir,'/voxel_map_0.20.mat'), 'file')==2
    map_gt = getfield(load(strcat(map_dir,'/voxel_map_0.20.mat')), 'map_gt');
else
    map_gt = VoxelMap(voxel_size, free_update, hit_update, occupied_threshold);
    map_gt = build_voxel_map(map_gt, map_dir, velo_to_map.T, measurement_range);
    save(strcat(map_dir,'/voxel_map_0.20.mat'), 'map_gt')
end

% Init maps: measurement, confidence, merged, ...
map_meas = VoxelMap(voxel_size, free_update, hit_update, occupied_threshold);
map_conf = VoxelMap(voxel_size, free_update, hit_update, occupied_threshold);
map_merged = VoxelMap(voxel_size, -100, 100, occupied_threshold);
map_measurable = VoxelMap(voxel_size, free_update, hit_update, occupied_threshold);

% Create voxelgrid in velo frame.
map_size = [320 320 32];
points_velo = gen_velo_points(map_size, -10);

% Create ray directions in sensor frame.
h_fov = [-pi/3 pi/3];
v_fov = [-pi/4 pi/4];
lidar_res = [160 120];
local_dirs = fov_dirs(h_fov, v_fov, lidar_res);

net = dagnn.DagNN.loadobj(getfield(load(net_path), 'net'));
net.conserveMemory = false;
if gpuDeviceCount() > 0
    net.move('gpu');
end

% Number of selected rays per position.
n_sel = 200;
% Position step.
step = 5;
% Planning horizont
n_plan_pos = 5;

if draw
    f1 = figure('Name', 'Measurements');
    grid on; axis equal;
    f2 = figure('Name', 'Occupancy Reconstruction');
    grid on; axis equal;
    f3 = figure('Name', 'Ground Truth');
    grid on; axis equal;
end

for frame = 1:step:size(path, 2)
    t_frame = tic();
    % Plan (random / greedy). Start with greedy.
    if strcmpi(planning, 'random') || frame < 1 + step
        i_rays = plan_rays('random', n_sel, local_dirs);
    else
        plan_pos = frame:step:frame+step*(n_plan_pos-1);
        plan_pos = plan_pos(plan_pos <= numel(velo_to_map.T));
        i_rays = plan_rays(planning, n_sel, local_dirs, velo_to_map.T(plan_pos), map_merged);
    end
    
    % Measure: trace rays in ground-truth map.
    origin = path(:, frame);
    dirs_in_map = velo_to_map.T{frame}(1:3,1:3) * local_dirs;
    [pos, val] = map_gt.trace_rays(origin, dirs_in_map(:, i_rays), measurement_range(2));
    % Update maps with valid measurements.
    map_meas.update_lines(origin, pos(:, ~isnan(val)));
    map_merged.update_lines(origin, pos(:, ~isnan(val)));
    
    % Estimate local occupancy from current measurements using CNN.
    points_in_map = p2e(velo_to_map.T{frame} * points_velo);
    [~, val] = map_meas.get_voxels(points_in_map);
    % Change inputs to [-1 0 1].
    val = sign(val);
    val(isnan(val)) = 0;
    input = single(reshape(val, map_size));
    if gpuDeviceCount() > 0
        input = gpuArray(input);
    end
    net.eval({'input',input})
    output = gather(net.vars(net.getVarIndex('output')).value);
    
    % Update global maps with local estimates.
    map_conf.update_voxels(points_in_map, double(output(:)'));
    map_merged.update_voxels(points_in_map, double(output(:)'));
    
    % Update measurable voxels.
    [pos,val] = map_gt.trace_rays(origin, dirs_in_map, measurement_range(2));
    map_measurable.update_lines(origin, pos(:, val >= occupied_threshold));
    
    fprintf('Frame %i / %i: %.3f s.\n', frame, size(path, 2), toc(t_frame));
    
    if draw
        clf(f1);
        plot_map(f1, map_meas, path(:,1:step:frame));
        drawnow;
        clf(f2);
        plot_map(f2, map_conf, path(:,1:step:frame));
        drawnow;
        clf(f3);
        plot_map(f3, map_gt, path(:,1:step:frame));
        drawnow;
    end
end

function [map_merge, map_conf, map_meas, map_gt, map_pos_input, path] = run_pipeline(map_name, net_name, planning, draw,mode)
% run matconvnet/matlab/vl_setupnn.m ; % activate MatConvNet if needed


%% init
base_dir = strcat('/datagrid/vras/petrito1/workspace/data/kitti/',map_name);
net_path = strcat('',net_name);

%% voxel map parameters
voxel_size = 0.2;
free_update = -1;
hit_update = 2;
occupied_threshold = 2;
measurement_range = [2.5 100];

%%
map_merge = VoxelMap(voxel_size, -100, 100, occupied_threshold);


% Load transforms from velo to map, create sensor path.
velo_to_map.T = load_velo_to_map(base_dir);
path = cell2mat(cellfun(@(T) T(1:3,4), velo_to_map.T, 'UniformOutput', false));

% build GT map
if exist(strcat(base_dir,'/voxel_map_0.20.mat'), 'file')==2
    map_gt = getfield(load(strcat(base_dir,'/voxel_map_0.20.mat')), 'vox_map');
else
    map_gt = VoxelMap(voxel_size, free_update, hit_update, occupied_threshold);
    map_gt = build_voxel_map(map_gt, base_dir, velo_to_map.T, measurement_range);
end

% Init maps: confidence, measurement, ...
map_conf = VoxelMap(voxel_size, free_update, hit_update, occupied_threshold);
map_meas = VoxelMap(voxel_size, free_update, hit_update, occupied_threshold);
map_pos_input = VoxelMap(voxel_size, free_update, hit_update, occupied_threshold);

% Create voxelgrid in velo frame.
map_size = [320 320 32];
points_velo = gen_velo_points(map_size, -10);

% Create ray directions in sensor frame.
h_fov = [-pi/3 pi/3];
v_fov = [-pi/4 pi/4];
lidar_res = [160 120];
local_dir = fov_dirs(h_fov, v_fov, lidar_res);

net = dagnn.DagNN.loadobj(getfield(load(net_path), 'net'));
net.conserveMemory = false;
if gpuDeviceCount() > 0
    net.move('gpu');
end

% Number of all possible rays per position.
n_rays = size(local_dir, 2);
% Number of selected rays per position.
n_sel = 200;
% Position step.
step = 5;
% Planning horizont
n_plan_pos = 5;

mkdir('map_exp');
if(draw)    
    f1 = figure(101);
    f2 = figure(102);
    f3 = figure(103);
end

%% pipeline
for frame=1:step:size(path, 2)
    t_frame = tic();
    %% Select measurement directions randomly or via planning.
    % Always start with random plan to get initial measurements.
    if strcmpi(planning, 'random') || frame < 1 + step
        i_rays = plan_rays('random', n_sel, local_dir);
    else
        plan_pos = frame:step:frame+step*(n_plan_pos-1);
        plan_pos = plan_pos(plan_pos <= numel(velo_to_map.T));
        i_rays = plan_rays(planning, n_sel, local_dir, velo_to_map.T(plan_pos), map_merge, mode);
    end
    
    direction = local_dir(:, i_rays);
    direction = velo_to_map.T{frame}(1:3,1:3)*direction;
    origin = path(:,frame);
    
    % trace rays in GT map
    [pos, val] = map_gt.trace_rays(origin, direction, measurement_range(2));
    
    points_in_map = p2e(velo_to_map.T{frame}*points_velo);
    % add input to M_meas
    map_meas.update_lines(origin, pos(:, ~isnan(val)));
    map_merge.update_lines(origin, pos(:, ~isnan(val)));
    
    %% merge
    % get CNN input
    [~, val] = map_meas.get_voxels(points_in_map);
    % upravit vstupy na -1 1 0
    val(isnan(val)) = 0;
    val(val<0) = -1;
    val(val>0) = 1;
    input = single(reshape(val, map_size));
    
    if gpuDeviceCount() > 0
        input = gpuArray(input);
    end
    
    % get CNN prediction
    net.eval({'input',input})
    output = gather(net.vars(net.getVarIndex('output')).value);
    
    map_conf.update_voxels(points_in_map, double(output(:)'));
    map_merge.update_voxels(points_in_map, double(output(:)'));
    
    %   measurable data
    direction = velo_to_map.T{frame}(1:3,1:3)*local_dir;
    [pos,val] = map_gt.trace_rays(origin, direction, measurement_range(2));
    map_pos_input.update_lines(origin, pos(:, val >= occupied_threshold));
    
    fprintf('Frame %i / %i: %.3f s.\n', frame, size(path, 2), toc(t_frame));
    
    if(draw==1)
        %% plot prediciton
        clf(f1)
        clf(f2)
        clf(f3)
        
        plot_map(f1, map_meas, path(:,1:step:frame));
        grid on
        drawnow
        
        plot_map(f2, map_conf, path(:,1:step:frame));
        grid on
        drawnow
        
        plot_map(f3, map_gt, path(:,1:step:frame));
        grid on
        drawnow
        axis equal
        
    end
    
    
end
function [A, c] = create_planning_task(map, T, local_dirs, mode)
%CREATE_PLANNING_TASK
%
% [A, c] = create_planning_task(map, T, local_dir, prob_vis)
%
% Input:
% - map Voxel map with logodds occupancy probability,
%       fused with the measurements so far.
% - T 1-by-P cell array with 4-by-4 rigid transforms from sensor to map.
% - local_dirs 3-by-D Ray directions in the sensor frame.
% - mode Visibility mode from {'occ_thresh', 'vis_thresh', 'vis_prob'}.
%
% Output:
% - A M-by-N sparse visibility matrix, where N = P*D.
% - c 1-by-M voxel values (gains).
%
assert(isa(map, 'VoxelMap'));
assert(iscell(T));
assert(ismatrix(local_dirs));
if nargin < 4 || isempty(mode)
    mode = 'vis_prob';
end
assert(ischar(mode));
assert(ismember(mode, {'occ_thresh' 'free_prob' 'vis_thresh' 'vis_prob' 'vis_prob_simple'}));
% mode = 'occ_thresh';  % A(i) from {0, 1}, all preceding voxels have occupancy within limits.
% mode = 'vis_thresh';  % A(i) from {0, 1}, all preceding voxels have visibility within limits.
% mode = 'vis_prob';    % A(i) from [0, 1].

n_pos = numel(T);
n_rays = size(local_dirs, 2);
prior = 1/12;  % Occupancy prior (approx.).
% thresh_logodds = logit(0.95);
thresh_logodds = logit(0.5);

%% Set ray tracing parameters.
% For map merged from estimates, all voxels within local neighborhood should be known.
% Otherwise min/max values should be nans and range should be reasonably low.
min_val = -inf;
if strcmpi(mode, 'occ_thresh')
    max_val = thresh_logodds;
else
    % max_val = inf;
    % prior + post < logit(0.999)
    max_val = logit(0.999) - logit(prior);
end
max_range = 48;
% max_range = inf;

velo_path = cell2mat(cellfun(@(T) T(1:3,4), T, 'UniformOutput', false));

%% Accumulate ray voxels for all position-directions.
t = tic();
ray_x = cell(1, n_pos);
ray_val = cell(1, n_pos);
vis_x = cell(1, n_pos);
for i_pos = 1:n_pos
    dirs = T{i_pos}(1:3, 1:3) * local_dirs;
    [~, ~, ray_x{i_pos}, ray_val{i_pos}] = ...
        map.trace_rays(velo_path(:, i_pos), dirs, max_range, min_val, max_val);
    % Get all visible voxels, per position, concatenate later.
    vis_x{i_pos} = cell2mat(ray_x{i_pos});
end
% fprintf('Ray tracing for %i positions: %.3f s.\n', n_pos, toc(t));
clear dirs from pos t;

%% Convert ray voxels to correseponing linear indices.
% Use row indices of sorted unique voxels for linear indices.
t = tic();
uni_x = unique(cell2mat(vis_x)', 'rows');
n_vox = size(uni_x, 1);
% Get voxel values for these voxels (logodds of occ. prob.).
[~, l] = map.get_voxels(uni_x');
% End points may be unknown/nan.
l = logit(prior) + l;
l(isnan(l)) = logit(prior);
assert(~any(isnan(l)));
p = logistic(l);
% c = gain(p);
% c = logit_entropy(l);
c = binary_entropy(p);
% fprintf('Gathering %i unique voxels: %.3f s.\n', n_vox, toc(t));

t = tic();
vox_ind = HashMap();
vox_ind.set(uni_x, (1:n_vox)');
% Rows and cols for sparse visibility matrix.
% Rows correspond to voxels, columns to position-rays.
rows_cell = cell(n_pos*n_rays, 1);
cols_cell = cell(n_pos*n_rays, 1);
% Index of a planned position, real ones can be skipped.
for i_pos = 1:n_pos
    for i_ray = 1:n_rays
        i_vox = vox_ind.get(ray_x{i_pos}{i_ray}');
        i_pos_ray = (i_pos-1)*n_rays + i_ray;
        rows_cell{i_pos_ray} = i_vox;
        cols_cell{i_pos_ray} = repmat(i_pos_ray, size(i_vox));
    end
end
% fprintf('Gathering visibility indices for %i positions: %.3f s.\n', n_pos, toc(t));
clear i_pos i_pos_ray i_ray i_vox t;

t = tic();
if strcmpi(mode, 'occ_thresh')
    %% Create and return tresholded visibility task instance.
    rows = cell2mat(rows_cell);
    cols = cell2mat(cols_cell);
    A = sparse(rows, cols, ones(size(rows)), n_vox, n_pos * n_rays);
    assert(all(sum(A, 2) > 0));
elseif strcmpi(mode, 'free_prob')
    rows = cell2mat(rows_cell);
    cols = cell2mat(cols_cell);
    % assert(all(size(p) == size(rows)));
    A = sparse(rows, cols, 1 - p(rows), n_vox, n_pos * n_rays);
    if ~all(sum(A, 2) > 0)
        fprintf('Zero row sums of A (visibility): %f.\n', mean(full(sum(A, 2)) == 0));
    end
    % assert(all(sum(A, 2) > 0));
else
    %% Create and return probabilistic visibility task instance.
    simple = strcmpi(mode, 'vis_prob_simple');
    P = prob_visibility(p, rows_cell, simple);
    % Prob. of current voxel being visible >= prob. of the next being
    % visible, so the order doesn't matter and we can safely threshold.
    if strcmpi(mode, 'vis_thresh')
        A = double(P >= logistic(thresh_logodds));
    else
        A = P;
    end
end

end

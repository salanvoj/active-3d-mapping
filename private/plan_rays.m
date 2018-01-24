function i_rays = plan_rays(plan, n_sel, dirs, T, map)
%PLAN_RAYS
%
% i_rays = plan_rays(plan, n_sel, local_dir, T,  map)
%
if isscalar(dirs)
    n_rays = dirs;
else
    n_rays = size(dirs, 2);
end

if strcmpi(plan, 'random')
    i_rays = randperm(n_rays, n_sel);
else
    t_task = tic();
    n_pos = numel(T);
    [A, c] = create_planning_task(map, T, dirs);
    t_task = toc(t_task);
    t_plan = tic();
    if strcmpi(plan, 'greedy')
        [x_ind, val, vox] = plan_rays_greedy(A, c, n_pos, n_sel);
        valid = x_ind < n_pos * n_rays;
        x_ind = x_ind(valid);
        x = zeros(n_pos * n_rays, 1);
    else
        error('Invalid planning type: %s', plan);
    end
    x(x_ind + 1) = 1;  % 0-based to 1-based indices.
    t_plan = toc(t_plan);
    vox = (1 - prod(1 - full(A(:, x>0)), 2))';
    val = sum(vox .* c);
    fprintf('Plan %s: %i rays, gain %.1f / %.1f, vox. %.1f / %i, %.3f s (task %.3f s).\n', ...
        plan, sum(x(1:n_rays)), val, sum(c(:)), sum(vox), numel(vox), t_plan, t_task);
    x = fill_plan(x, n_pos, n_rays, n_sel);
    i_rays = find(x(1:n_rays));
    assert(numel(i_rays) == n_sel);
end

end

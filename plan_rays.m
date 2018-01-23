function i_rays = plan_rays(plan, n_sel, dirs, T, map, mode)
%PLAN_RAYS Summary of this function goes here
%
% i_rays = plan_rays(plan, n_sel, local_dir, T,  map, mode)
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
    [A, c] = create_planning_task(map, T, dirs, mode);
%     fprintf('Creating planning task (%i nnz): %.3f s.\n', nnz(A), toc(t));
    t_task = toc(t_task);
    t_plan = tic();
    if strcmpi(plan, 'greedy')
%         [x, val, vox] = plan_rays_greedy(A, c, n_pos, n_sel, false);
        [x, val, vox] = plan_rays_greedy(A, c, n_pos, n_sel, true);
%         [x, val, vox] = plan_rays_greedy(A, c, n_pos, n_sel, 2);
        % val and vox need to be recomputed.
    elseif strcmpi(plan, 'lp')
        [x, val, vox] = plan_rays_lp(A, c, n_pos, n_sel);
    elseif strcmpi(plan, 'milp')
        [x, val, vox] = plan_rays_milp(A, c, n_pos, n_sel);
    end
    t_plan = toc(t_plan);
%     fid = fopen('/tmp/active_mapping/t_task.txt', 'a');
%     fprintf(fid, '%.6f\n', t_task);
%     fclose(fid);
%     fid = fopen('/tmp/active_mapping/t_plan.txt', 'a');
%     fprintf(fid, '%.6f\n', t_plan);
%     fclose(fid);
    vox = (1 - prod(1 - full(A(:, x>0)), 2))';
    val = sum(vox .* c);
    fprintf('Plan %s: %i rays, gain %.1f / %.1f, vox. %.1f / %i, %.3f s (task %.3f s).\n', ...
        plan, sum(x(1:n_rays)), val, sum(c(:)), sum(vox), numel(vox), t_plan, t_task);
    % check_plan(planning, x_sol, n_pos, n_sel, A, c);
    x = fill_plan(x, n_pos, n_rays, n_sel);
    i_rays = find(x(1:n_rays));
    assert(numel(i_rays) == n_sel);
end

end

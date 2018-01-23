function [x, val, vox] = plan_rays_greedy(A, c, n_pos, n_sel)
%PLAN_RAYS_GREEDY Plan rays for maximum weighted coverage greedy ray by ray
%
% [x, val, vox] = plan_rays_greedy(Vis, vox_val, n_pos, n_sel)
%
% Input:
% - A M-by-N Visibility matrix where M is the number of voxels and N = P*Q where
%   P is the number of positions and Q the number of viewpoints.
% - c 1-by-M or M-by-1 voxel value (gain) vector.
% - n_pos Number of positions.
% - n_sel Number of selected views per position.
%
% Output:
% - x N-by-1 vector of selected rays for all positions.
% - val Weighted voxel coverage.
% - vox 1-by-M voxel coverage mask, covered = true, not covered = false.

% Check types (coder) and consistency.
% assert(issparse(A) && isa(A, 'double') && size(A, 1) >= 1 && size(A, 2) >= 1);
% assert(size(c, 1) == 1 && size(c, 2) >= 1);
assert(ismatrix(A) && issparse(A));
assert(all(A(A ~= 0) == 1));
% Use sparse double matrix, to be able to represent ray/vox visibility probability.
A = double(A);
assert(isvector(c) && isa(c, 'double'));
assert(size(A, 1) == numel(c));
assert(all(c >= 0));
c = c(:);
assert(isscalar(n_pos) && isa(n_pos, 'double'));
assert(isscalar(n_sel) && isa(n_sel, 'double'));

t0 = tic();
t = t0;

n_rays = size(A, 2) / n_pos;
sel = zeros([1 n_pos]);
x = false([n_pos*n_rays 1]);
% Ac = bsxfun(@times, A, c);  % too slow
[i_, j_, v_] = find(A);
C = sparse(i_, j_, v_ .* c(i_), size(A, 1), size(A, 2));
clear i_ j_ v_;
while any(sel < n_sel)
    [~, i] = max(sum(C));
    C(:, i) = 0;
    pos = floor((i - 1) / n_rays) + 1;
    assert(pos >= 1 && pos <= n_pos);
    if sel(pos) == n_sel
        % All voxels are already covered.
        % Missing rays will be selected in fix_plan.
        break;
    end
    % If there are still rays to be selected at the given position.
    if sel(pos) < n_sel
        x(i) = true;
        % Zero contributions for the visible voxels for all rays.
        C(A(:, i) > 0, :) = 0;  % maybe slightly slower
%         C(find(A(:, i)), :) = 0;
        sel(pos) = sel(pos) + 1;
        if sel(pos) == n_sel
            % Zero out all rays for given position.
            C(:, (pos-1)*n_rays+1:pos*n_rays) = 0;
        end
    end
    if toc(t) >= 30
        fprintf('Greedy planning: %.1f %% (%.1f s).\n', 100 * mean(sel) / n_sel, toc(t0));
        t = tic();
    end
end

% x = fix_plan(x, n_pos, n_rays, n_sel);
val = c' * (A * x > 0);
vox = A * x > 0;

end

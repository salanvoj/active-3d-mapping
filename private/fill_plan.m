function x = fill_plan(x, n_pos, n_rays, n_sel)
%FILL_PLAN Select unused measurements randomly
%
% x = fill_plan(x, n_pos, n_rays, n_sel)
%
assert(isvector(x));
assert(isscalar(n_pos));
assert(isscalar(n_rays));
assert(isscalar(n_sel));
assert(n_sel <= n_rays);
assert(numel(x) >= n_pos * n_rays);

x = logical(x(:));
x = x(1:n_pos*n_rays);
for i = 1:n_pos
    % Get ray indices for position i.
    i_rays = (i-1)*n_rays+1:i*n_rays;
    n_fill = n_sel - sum(x(i_rays));
    % assert(n_fill >= 0);
    if n_fill < 0
        warning('Position %i, planned rays: %.1f from %i.', i, sum(x(i_rays)), n_sel);
        n_fill = 0;
    end
    if n_fill == 0
        continue;
    end
    i_fill = i_rays(~x(i_rays));
    i_fill = i_fill(randperm(numel(i_fill), n_fill));
    assert(numel(i_fill) == n_fill);
    x(i_fill) = true;
end
x = sparse(x);
assert(nnz(x) == n_pos * n_sel);

end

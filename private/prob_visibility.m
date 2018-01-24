function A = prob_visibility(p, r, simple)
%PROB_VISIBILITY Probabilistic ray-voxel visibility matrix
%
% A = prob_visibility(p, r)
%
% Input:
% - p M-by-1 Voxel occupancy probability, for M voxels.
% - r N-by-1 cell array of voxel indices, for N rays.
% - simple logical Need an occupied voxel behind?
%
% Output:
% - A M-by-N probabilistic visibility matrix, where
%   A(i,j) is the probability of visibility of the ith voxel
%   from the jth ray.
%

if nargin < 3 || isempty(simple)
    simple = false;
end

assert(all(p >= 0));
assert(all(p <= 1));
assert(iscell(r));

p = p(:);
m = numel(p);
n = numel(r);

% Arrange everything in columns.
r = r(:);

c = cell(size(r));
v = cell(size(r));
for i = 1:numel(r)
    r{i} = r{i}(:);
    c{i} = repmat(i, size(r{i}));
    % Prob. of the voxel being visible/measurable is
    % the prob. of all the preceding voxels being empty times
    % the prob. of any voxel following within given range being occupied.
    vf = cumprod(1 - p(r{i}));
    if simple
        v{i} = [1; vf((1:end-1)')];
        % fprintf('Simple visibility prob. computation.\n');
    else
        vr = cumprod(1 - p(r{i}), 'reverse');
        v{i} = [1; vf((1:end-1)')] .* (1 - vr);
    end
end
r = cell2mat(r);
c = cell2mat(c);
v = cell2mat(v);

A = sparse(r, c, v, m, n);

end

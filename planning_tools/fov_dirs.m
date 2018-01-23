function dirs = fov_dirs(h_fov, v_fov, res)
%FOV_DIRS Uniformly distributed directions in given field of view
%
% dirs = fov_dirs(h_fov, v_fov, res)
%
% - h_fov  Horizontal field of view, scalar angle or interval [min max],
%          in radians.
% - v_fov  Vertical field of view, scalar angle or interval [min max],
%          in radians.
% - res    Sensor resolution, [width height].
% 
% - dirs  3-by-prod(res) matrix of ray directions.
%
if isscalar(h_fov)
    h_fov = [-h_fov/2 h_fov/2];
end
if isscalar(v_fov)
    v_fov = [-v_fov/2 v_fov/2];
end
if isscalar(res)
    res = [res res];
end

[az, el] = ndgrid(linspace(h_fov(1), h_fov(2), res(1)), ...
                  linspace(v_fov(1), v_fov(2), res(2)));
[dx, dy, dz] = sph2cart(az, el, 1);
clear az el;
dirs = [dx(:)'; dy(:)'; dz(:)'];

end

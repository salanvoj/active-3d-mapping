function points_velo = gen_velo_points(mat_size,x_shift)
% generates local map (grid) in velodine frame
% mat_size - matrix size
% x_shift - shift of velodine sensor in x direction


[x, y, z] = ndgrid(linspace(-31.9-(x_shift), 31.9-(x_shift), mat_size(1)), ...
    linspace(-31.9, 31.9, mat_size(2)), ...
    linspace( 3.6, -2.6, mat_size(3)));

points_velo = e2p([x(:) y(:) z(:)]');

end
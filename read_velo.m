function [x, i] = read_velo(path, range)
%READ_VELO Read velodyne points
%
% [x, i] = read_velo(cloud_path)
% [x, i] = read_velo(cloud_path, range)
%
if nargin < 2 || isempty(range)
    range = [0 inf];
end
assert(isa(range, 'double'));
assert(numel(range) == 2);

fid = fopen(path, 'rb');
x = fread(fid, [4 inf], 'single')';  % In velodyne coordinates
fclose(fid);
i = x(:, 4)';    % Intensity
x = x(:, 1:3)';  % Points in columns

% Filter out points to close to the sensor (on the car)
% and points which are too far.
x_norm = sqrt(sum(x.^2));
keep = x_norm >= range(1) & x_norm <= range(2);
x = x(:, keep);
end

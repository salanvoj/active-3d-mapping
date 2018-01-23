function T = load_velo_to_map(base_dir)
%LOAD_VELO_TO_MAP Load velo to map transforms
%
% T = load_velo_to_map(base_dir)
%
% Input:
% - base_dir Directory containing oxts.
%
% Output:
% - T 1-by-N cell array with 4-by-4 rigid transforms.
%

oxts = loadOxtsliteData(base_dir);
T = convertOxtsToPose(oxts);
T_imu_to_velo = loadCalibrationRigid(fullfile(base_dir, '..', 'calib_imu_to_velo.txt'));
T = cellfun(@(T) T / T_imu_to_velo, T, 'UniformOutput', false);

end
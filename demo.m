mex_planning();  % Comment out once compiled.
close all; clear; clc;

% testing sequences
map_name{1} = '2011_09_26/2011_09_26_drive_0048_sync';
map_name{2} = '2011_09_26/2011_09_26_drive_0095_sync';
map_name{3} = '2011_09_26/2011_09_26_drive_0113_sync';
map_name{4} = '2011_09_26/2011_09_26_drive_0001_sync';
map_name{5} = '2011_09_26/2011_09_26_drive_0002_sync';
map_name{6} = '2011_09_26/2011_09_26_drive_0005_sync';
map_name{7} = '2011_09_26/2011_09_26_drive_0060_sync';
map_name{8} = '2011_09_26/2011_09_26_drive_0091_sync';
map_name{9} = '2011_09_26/2011_09_26_drive_0009_sync';
map_name{10} = '2011_09_26/2011_09_26_drive_0084_sync';
map_name{11} = '2011_09_26/2011_09_26_drive_0095_sync';
map_name{12} = '2011_09_26/2011_09_26_drive_0104_sync';
map_name{13} = '2011_09_26/2011_09_26_drive_0106_sync';
% validation sequences
map_name{14} = '2011_09_26/2011_09_26_drive_0039_sync';
map_name{15} = '2011_09_26/2011_09_26_drive_0046_sync';
map_name{16} = '2011_09_26/2011_09_26_drive_0020_sync';
% training sequences
map_name{17} = '2011_09_26/2011_09_26_drive_0019_sync';
map_name{18} = '2011_09_26/2011_09_26_drive_0022_sync';
map_name{19} = '2011_09_26/2011_09_26_drive_0023_sync';
map_name{20} = '2011_09_26/2011_09_26_drive_0035_sync';
map_name{21} = '2011_09_26/2011_09_26_drive_0036_sync';
map_name{22} = '2011_09_26/2011_09_26_drive_0061_sync';
map_name{23} = '2011_09_26/2011_09_26_drive_0064_sync';
map_name{24} = '2011_09_26/2011_09_26_drive_0079_sync';
map_name{25} = '2011_09_26/2011_09_26_drive_0086_sync';
map_name{26} = '2011_09_26/2011_09_26_drive_0087_sync';

index = 1; % change index to choose map;  
choosen_map = map_name{index};

net_name = 'final_net.mat';
planning = 'greedy';  % change to 'random' to make random measurements
draw_figure = true;  % change to false to disable drawing figures

[map_merged, map_conf, map_meas, map_gt, map_measurable, path] = ...
    run_pipeline(choosen_map, net_name, planning, draw_figure);

clc;
clear all;
close all;

%% Graph and Find Flags
[flag_1, time_tle_1, r_output_1, v_output_1] = Flag_Finder(31136);
%[flag_2, time_tle_2, r_output_2, v_output_2] = Flag_Finder_Updated();

% 1. debug Flag_Finder_updated and Orbit_Prop_Test (seems to work but need to
    % waitn long enough for plots to show)

% 2. find the union times of two different satellites (the times that they
% are both interpolated to, the indexes of the r positions at those times
% for eah satellite, and create the filtered r vectors for each satellite

% 3. Use eciToLvlhRelativeStates.m function to get the relative distances
% at different times

% 4.  Plot the relative dsitancaes

%% 
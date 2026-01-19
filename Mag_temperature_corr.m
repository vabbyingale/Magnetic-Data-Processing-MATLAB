function [XYZ_tcal,mag_tcal] = Mag_temperature_corr(sensor_number, mag_x, mag_y, mag_z, mag_temp, nbr_ind)
    % A function for temperature correction to magnetometer data using CALIB1520 standard.
    %
    % Syntax:
    %   [XYZ_tcal, mag_tcal] = temperature_corr(690, mag.x, mag.y, mag.z, mag.temp, find(idx));
    %
    % Inputs:
    %   sensor_number : Integer specifying sensor ID (e.g., 690 for top fluxgate)
    %   mag_x, mag_y, mag_z : Vectors of raw magnetic field measurements
    %   mag_temp      : Vector of temperature measurements
    %   nbr_ind       : Indices of samples to process
    %
    % Outputs:
    %   XYZ_tcal      : Nx3 matrix of temperature-corrected magnetic data [X Y Z]
    %   mag_tcal      : Full output from calib1520_sensors (may include all corrections)


% Define variables of magnetic data for duration of segment
mag_xt    = mag_x(nbr_ind)';
mag_yt    = mag_y(nbr_ind)';
mag_zt    = mag_z(nbr_ind)';
mag_tempt = mag_temp(nbr_ind)';
mag_xyzt  = [mag_xt,mag_yt,mag_zt,mag_tempt];

% Sesnor information is from calib1520_sensors code 
% sensor number is 0690 for top fluxgate magnetometer

[XYZ_tcal,mag_tcal] = Mag_calib1520_sensors(mag_xyzt,sensor_number);


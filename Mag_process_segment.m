function Mag_process_segment(mag, sbe, Fitting)
    % A function to process a selected straight-line segment of magnetic data.
    % Prompts user for a segment duration, finds the straight segment.
    %
    % Syntax:
    %   process_segment(mag, sbe, Fitting)
    %
    % Inputs:
    %   mag     : Struct with magnetometer data (utime, x, y, z, heading, pitch, roll, temp, lat, lon, samp_rate)
    %   sbe     : Struct with SBE data (utime, lat, lon) used for segment finding
    %   Fitting : Instrument-specific heading correction parameters
    % 
    % Calls out:
    %   find_straight_segment.m, plot_filter_components.m, magfd.m
    
close all;
%-------------------------------------------------------------------------
% Step 1: user defined segment duration

duration = input('Segment duration (s): ');
[seg_start, seg_end] = Mag_find_straight_segment(duration, sbe.utime, sbe.lon, sbe.lat, mag.utime);
idx = mag.utime >= seg_start & mag.utime <= seg_end;

%-------------------------------------------------------------------------
% Step 2: extract variables
mag_ut = mag.utime(idx); mag_dt = datetime(mag_ut,'ConvertFrom','epochtime','Epoch','01-Jan-1970');
mag_x = mag.x(idx)*1e5; mag_y = mag.y(idx)*1e5; mag_z = mag.z(idx)*1e5;
hdg = mag.heading(idx); pitch = mag.pitch(idx); roll = mag.roll(idx);

%-------------------------------------------------------------------------
% Step 3: Analysis

% Plotting
win = round(5*(mag.samp_rate/0.0166667));       % For filtering
Mag_plot_filter_components(mag_dt, mag_x, mag_y, mag_z, win,'Mag X','Mag Y','Mag Z',1);

% Temperature correction
sensor = 690;
[XYZ_tcal,~] = Mag_temperature_corr(sensor, mag.x, mag.y, mag.z, mag.temp, find(idx));

% Component rotation
[ned] = Mag_component_rotation(hdg,pitch,roll,XYZ_tcal,0);
n = ned(:,1); e = ned(:,2); d = ned(:,3);

% Plot filter components
Mag_plot_filter_components(mag_dt, n,e,d, win,'North','East','Down',1);
corn = Fitting(1)+Fitting(4)*sind(hdg+Fitting(7));
core = Fitting(2)+Fitting(5)*sind(hdg+Fitting(8));
corz = Fitting(3)+Fitting(6)*sind(hdg+Fitting(9));

% Heading correction
Nc = n - corn; Ec = e - core; Dc = d - corz;
Mag_corr = Mag_plot_filter_components(mag_dt, Nc,Ec,Dc, win,'North-corrected','East-corrected','Down-corrected',2);
mag_total = sqrt(sum(Mag_corr.^2,2));

% IGRF model
lat1 = mag.lat(find(idx, 1));  % Get first index where idx is true
lon1 = mag.lon(find(idx, 1));
J = Mag_magfd(2024.9, 1, 0.003, 90 - lat1, lon1);  % 0.003 km ≈ 3 m altitude

IGRF = [J(1), J(2), J(3)];
Dec = 11.0707; Inc = -21.5330;

% Anomaly computation
anom_n = Nc - IGRF(1); anom_e = Ec - IGRF(2); anom_d = Dc - IGRF(3);
anom_t = anom_n * cosd(Dec) * cosd(Inc) + anom_e * sind(Dec) * cos(Inc) + anom_d * sin(Inc);
Mag_plot_filter_components(mag_dt, anom_n,anom_e,anom_d, win,'North-Anomaly','East-Anomaly','Down-Anomaly',0);

function Mag_process_dive_without_spins(mag, sbe, Fitting, file_prefix)
    % A funcion to process and save magnetic data exlcuding spins
    % Applies spin correction, heading corrections, filters, computes anomaly,
    % saves to .mat and .txt files, with header.
    %
    % Syntax:
    %   process_segment(mag, sbe, Fitting, file_prefix)
    %
    % Inputs:
    %   mag         - Struct with magnetometer data
    %   sbe         - Struct with SBE data
    %   Fitting     - Spin fitting coefficients
    %   file_prefix - File name prefix for outputs
    %
    % Outputs:
    %   Saves .mat and .txt files with processed data
    %
    % Calls out:
    %   dive_segment_entire.m, temperature_corr.m, component_rotation.m,
    %   plot_filter_components.m. magfd.m

close all;
%-------------------------------------------------------------------------
% Step 1: Filename initiation and checking

matFile = sprintf('Seg_%s.mat', file_prefix);
txtFile = sprintf('Seg_%s.txt', file_prefix);

% Check if file already exists
if isfile(matFile)
    fprintf('\n')
    fprintf('File already exists: %s\n', matFile);
    choice = input('Recompute anyway? (y/n): ','s');
    if lower(choice) ~= 'y'
        fprintf('Check for already saved file.\n');
        return;
    else
        fprintf('Recomputing as requested.\n');
    end
end

%-------------------------------------------------------------------------
% Step 2: segment selection
[seg_start, seg_end] = Mag_select_dive_segment(sbe.heading, sbe.lon, sbe.lat, sbe.utime, sbe.alt);
idx = mag.utime >= seg_start & mag.utime <= seg_end;
N = sum(idx);

%-------------------------------------------------------------------------
% Step 3: compute variables
mag_ut = mag.utime(idx);
mag_dt = datetime(mag_ut,'ConvertFrom','epochtime','Epoch','01-Jan-1970');
hdg = mag.heading(idx);
pitch = mag.pitch(idx);
roll = mag.roll(idx);

% Magnetometer data
mag_x = mag.x(idx) * 1e5;
mag_y = mag.y(idx) * 1e5;
mag_z = mag.z(idx) * 1e5;

%-------------------------------------------------------------------------
% Step 4: data processing

% Temperature correction
[XYZ_tcal, ~] = Mag_temperature_corr(690, mag.x, mag.y, mag.z, mag.temp, find(idx));

% Coordinate rotation
ned = Mag_component_rotation(hdg, pitch, roll, XYZ_tcal, 0);
n = ned(:,1);
e = ned(:,2);
d = ned(:,3);

% Heading-dependent correction
corn = Fitting(1) + Fitting(4) * sind(hdg + Fitting(7));
core = Fitting(2) + Fitting(5) * sind(hdg + Fitting(8));
corz = Fitting(3) + Fitting(6) * sind(hdg + Fitting(9));

Nc = n - corn;
Ec = e - core;
Dc = d - corz;

% Filtering
win = 3*round(mag.samp_rate / 0.0166667);
Filtered = Mag_plot_filter_components(mag_dt, Nc, Ec, Dc, win, 'North-corrected', 'East-corrected', 'Down-corrected', 1);

% Field strengths
mag_total_corr = sqrt(sum(Filtered.^2,2));
mag_total_orig = sqrt(n.^2 + e.^2 + d.^2);

% IGRF model
lat1 = mag.lat(find(idx,1));
lon1 = mag.lon(find(idx,1));
J = Mag_magfd(2024.9, 1, 0.003, 90 - lat1, lon1);
% J = magfd(2021.33, 1, 0.003, 90 - lat1, lon1);

% Declination & inclination
Dec = 11.0707;
Inc = -21.5330;

% Dec = 6.3639; Inc = 32.3010;

% Magnetic anomaly
anom = (Filtered(:,1)-J(1))*cosd(Dec)*cosd(Inc) + ...
       (Filtered(:,2)-J(2))*sind(Dec)*cosd(Inc) + ...
       (Filtered(:,3)-J(3))*sind(Inc);

% Resample SBE data
sbe_depth1 = resample(sbe.depth, mag.samp_rate, sbe.samp_rate);
sbe_temp1  = resample(sbe.temp,  mag.samp_rate, sbe.samp_rate);
sbe_alt1   = resample(sbe.alt,   mag.samp_rate, sbe.samp_rate);

% Truncate to match
sbe_depth1 = sbe_depth1(idx);
sbe_temp1  = sbe_temp1(idx);
sbe_alt1   = sbe_alt1(idx);

%-------------------------------------------------------------------------
% Step 4: save the files
output = [mag_ut', mag.lat(idx), mag.lon(idx), ...
          sbe_depth1, sbe_temp1', sbe_alt1, ...
          hdg, roll, pitch, ...
          mag_x', mag_y', mag_z', ...
          n, e, d, ...
          mag_total_orig, ...
          Filtered(:,1), Filtered(:,2), Filtered(:,3), ...
          mag_total_corr, ...
          Filtered(:,1)-J(1), Filtered(:,2)-J(2), Filtered(:,3)-J(3), ...
          anom];

header = 'utime,lat,lon,Depth,Temp,Altitude,Heading,Roll,Pitch,mag_x,mag_y,mag_z,N,E,D,mag_total_orig,N_corr,E_corr,D_corr,mag_total_corr,Anom_n,Anom_e,Anom_d,Anom_total';

% Downsample
output = downsample(output, 6);

% Save MAT
mag_top = struct();
mag_top.utime = output(:,1);
mag_top.lat   = output(:,2);
mag_top.lon   = output(:,3);
mag_top.depth = output(:,4);
mag_top.temp  = output(:,5);
mag_top.alt   = output(:,6);
mag_top.heading = output(:,7);
mag_top.roll    = output(:,8);
mag_top.pitch   = output(:,9);
mag_top.mag_x   = output(:,10);
mag_top.mag_y   = output(:,11);
mag_top.mag_z   = output(:,12);
mag_top.n       = output(:,13);
mag_top.e       = output(:,14);
mag_top.d       = output(:,15);
mag_top.mag_total_orig = output(:,16);
mag_top.N_corr  = output(:,17);
mag_top.E_corr  = output(:,18);
mag_top.D_corr  = output(:,19);
mag_top.mag_total_corr = output(:,20);
mag_top.anom_n = output(:,21);
mag_top.anom_e = output(:,22);
mag_top.anom_d = output(:,23);
mag_top.anom_total = output(:,24);

save(matFile, 'mag_top');

% Save TXT
fid = fopen(txtFile, 'w');
fprintf(fid, '%s\n', header);
fclose(fid);
dlmwrite(txtFile, output, '-append', 'precision', '%.6f');

fprintf('Dive data saved to %s and \n', matFile);
fprintf('                   %s     \n', txtFile);



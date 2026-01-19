function [input_for_upcont] = Mag_distance_along_track(sentry_number, date_part, time_part)
    % Computes track distance and saves formatted file as input for upward continuation
    %
    % Syntax: 
    %   input_for_upcont = distance_along_track(sentry_number, date_part, time_part)
    %
    % Inputs:
    %   sentry_number = String of dive number (e.g., '744')
    %   date_part     = String of dive date in yyyymmdd format (e.g., '20241204')
    %   time_part     = String of dive start time in hhmm format (e.g., '1654')
    %
    % Outputs:
    %   input_for_guspi - Nx8 matrix with formatted data:
    %       [distance_km, lon, lat, fish_depth_km, anom_n, anom_e, anom_d, anom_t]
    %
    % Notes:
    %   Requires a file 'Seg_sentry<dive>_<date>_<time>.txt' with navigation and anomaly data.
    %   The output '<base>_guspi.txt' is saved in the same directory.

% -------------------------------------------------------------------------
% Step 1: Build file names
base_prefix     = sprintf('Seg_sentry%s_%s_%s', sentry_number, date_part, time_part);
input_file      = sprintf('%s_cross.txt', base_prefix);
temp_dist_file  = sprintf('%s_distance.txt', base_prefix);
output_file     = sprintf('%s_upcont.txt', base_prefix);

% -------------------------------------------------------------------------
% Step 2: Load file and extract lat/lon
Path = readmatrix(input_file);
Lat = Path(:,2);
Lon = Path(:,3);

% -------------------------------------------------------------------------
% Step 3: Compute cumulative distance using Haversine formula
n_points = length(Lat);
distance = zeros(n_points,1);

for i = 2:n_points
    distance(i) = distance(i-1) + Mag_haversine(Lat(i-1), Lon(i-1), Lat(i), Lon(i));
end

% -------------------------------------------------------------------------
% Step 4: Save temporary distance file
writetable(table(distance), temp_dist_file, 'Delimiter', '\t', ...
           'FileType', 'text', 'WriteVariableNames', false);

% -------------------------------------------------------------------------
% Step 5: Reload full dataset and distance
data1 = readmatrix(input_file);
data2 = readmatrix(temp_dist_file);
distance = data2(:,1);      % distance in km  

% Extract field of interest
lat = data1(:,2);           % latitude
lon = data1(:,3);           % longitude
fish = data1(:,4) / 1000;   % depth in km
anom_n = data1(:,25);       % north-anomaly
anom_e = data1(:,26);       % east-anomaly
anom_d = data1(:,27);       % down-anomaly
anom_t = data1(:,28);       % total-anomaly

valid_idx = ~any(isnan([lat, lon, fish, anom_n, anom_e, anom_d, anom_t, distance]), 2);

% Apply filter to all relevant variables
lat      = lat(valid_idx);
lon      = lon(valid_idx);
fish     = fish(valid_idx);
anom_n   = anom_n(valid_idx);
anom_e   = anom_e(valid_idx);
anom_d   = anom_d(valid_idx);
anom_t   = anom_t(valid_idx);
distance = distance(valid_idx);

% -------------------------------------------------------------------------
% Step 6: Merge and write formatted output
merged_data = [distance, lon, lat, fish, anom_n, anom_e, anom_d, anom_t];

fid = fopen(output_file, 'w');
fmt = '%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n';

for i = 1:size(merged_data,1)
    fprintf(fid, fmt, merged_data(i,:));
end
fclose(fid);

% -------------------------------------------------------------------------
% Step 7: Clean up and return

delete(temp_dist_file);
input_for_upcont = merged_data;

fprintf('Input file for upward continuation is ready:\n %s\n', output_file);
end






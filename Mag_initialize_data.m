function [mag, sbe] = Mag_initialize_data(renav_mag, renav_sbe49)
    % intializeData: Initializes structured data from raw renav inputs
    %
    % Syntax:
    %   [mag, sbe] = initializeData(renav_mag, renav_sbe49);
    %
    % Inputs:
    %   renav_mag   = struct containing maggie_top fields (from re-navigation)
    %   renav_sbe49 = struct containing sbe49 fields (from re-navigation)
    %
    % Outputs:
    %   sbe = struct with time, position, depth, heading, altitude, temp, sampling rate
    %   mag = struct with time, mag xyz, attitude, position, temp, sampling rate

%-------------------------------------------------------------------------
% Step 1: sbe data

fprintf('\n')
disp('To verify variables, refer to function Mag_initialize_data')
sbe.utime     = renav_sbe49.t;
sbe.dtime     = datetime(sbe.utime,'ConvertFrom','epochtime','Epoch','01-Jan-1970');
sbe.lat       = renav_sbe49.renavLat;
sbe.lon       = renav_sbe49.renavLon;
sbe.depth     = renav_sbe49.renavDepth;
sbe.heading   = renav_sbe49.renavHeading;
sbe.alt       = renav_sbe49.renavAltitude;
sbe.temp      = renav_sbe49.temperature;
sbe.samp_rate = round(1/((sbe.utime(end)-sbe.utime(1))/numel(sbe.utime)));

%-------------------------------------------------------------------------
% Step 2: mag_renav data
mag.utime     = renav_mag.maggie_top.t;
mag.dtime     = datetime(mag.utime,'ConvertFrom','epochtime','Epoch','01-Jan-1970');
mag.x         = renav_mag.maggie_top.x;
mag.y         = renav_mag.maggie_top.y;
mag.z         = renav_mag.maggie_top.z;
mag.temp      = renav_mag.maggie_top.temperature;
mag.heading   = renav_mag.maggie_top.renavHeading;
mag.roll      = renav_mag.maggie_top.renavRoll;
mag.pitch     = renav_mag.maggie_top.renavPitch;
mag.lat       = renav_mag.maggie_top.renavLat;
mag.lon       = renav_mag.maggie_top.renavLon;
mag.samp_rate = round(1/((mag.utime(end)-mag.utime(1))/numel(mag.utime)));


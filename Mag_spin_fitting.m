function Fitting = Mag_spin_fitting(mag, sbe, spinFitFile, spinTimeFile, sentry_number)
    % Performs spin fitting calibration for a magnetometer survey
    % - Identifies a spin maneuver (period where vehicle rotated on spot)
    % - Applies temperature correction and coordinate rotation
    % - Fits heading-dependent corrections to remove spin-induced distortions
    % - Saves fitted calibration parameters and spin period times
    %
    % Syntax:
    %   Fitting = spin_fitting(mag, sbe, 'spin_fit.mat', 'spin_time.mat');
    %
    % Inputs:
    %   mag          : Struct of magnetometer data (utime, x, y, z, heading, pitch, roll, temp, samp_rate)
    %   sbe          : Struct of navigation data (heading, lon, lat, utime, alt)
    %   spinFitFile  : Output filename to save fitted calibration parameters
    %   spinTimeFile : Output filename to save spin start/end time
    %
    % Output:
    %   Fitting      : Vector of fitted calibration coefficients (used to correct heading-dependent errors)
    %

[spin_start, spin_end] = Mag_find_spins(sbe.heading, sbe.lon, sbe.lat, sbe.utime, sbe.alt, sentry_number);

% ------------------------------------------------------------------------
% Step 1: initialize variables

idx = mag.utime >= spin_start & mag.utime <= spin_end;
mag_ut = mag.utime(idx);
mag_dt = datetime(mag_ut,'ConvertFrom','epochtime','Epoch','01-Jan-1970');
mag_x = mag.x(idx)*1e5;                                 % gauss to nT
mag_y = mag.y(idx)*1e5;                                 % gauss to nT
mag_z = mag.z(idx)*1e5;                                 % gauss to nT

hdg = mag.heading(idx);
pitch = mag.pitch(idx);
roll = mag.roll(idx);

% ------------------------------------------------------------------------
% Step 2: analysis

% Filtering: 0.0166667 is the number to get total samples per minute
win = round(0.1*(mag.samp_rate/0.0166667)); 
Mag_plot_filter_components(mag_dt, mag_x, mag_y, mag_z, win,'Mag X','Mag Y','Mag Z',1);

% Temperature correction
sensor = 690;
[XYZ_tcal,~] = Mag_temperature_corr(sensor, mag.x, mag.y, mag.z, mag.temp, find(idx));
[ned] = Mag_component_rotation(hdg,pitch,roll,XYZ_tcal,0);
n = ned(:,1); e = ned(:,2); d = ned(:,3);
Mag_plot_filter_components(mag_dt, n,e,d, win,'North','East','Down',1);

% Getting spin data
Fitting = Mag_fit_spin_comp(hdg,n,e,d, sentry_number);
fprintf('#--------------- Spin calibration complete ---------------#\n');

% ------------------------------------------------------------------------
% Step 3: save spin  data

save(spinFitFile, 'Fitting');
spin_time = [spin_start; spin_end];
save(spinTimeFile, 'spin_time');

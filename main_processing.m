% main_processing.m - Entry point for Sentry magnetic data processing
%
% Processes magnetic and SBE data from Sentry dives, performs spin calibration,
% and saves processed data in multiple forms (segment, full dive, no spins).
%
% Author: Vaibhav Vijay Ingale

clc; clear all; close all;
warning('off', 'all');

%%

fprintf('===========================================================\n');
fprintf('           Sentry Magnetic Data Processing Pipeline        \n');
fprintf('===========================================================\n');

%-------------------------------------------------------------------------
% Step 1: Get user inputs for Sentry dive identification

fprintf('\n======================== Load the file ====================\n');
sentry_number = input('Enter Sentry number (e.g. 744): ', 's');
date_part     = input('Enter date (yyyymmdd, e.g. 20241204): ', 's');  
time_part     = input('Enter file start time (e.g. 1654): ', 's');

% Build prefix used in filenames
file_prefix = sprintf('sentry%s_%s_%s', sentry_number, date_part, time_part);

% Input data files
sbe_file = sprintf('%s_sbe49_renav.mat', file_prefix);  % SBE sensor data
mag_file = sprintf('%s_mag_renav.mat', file_prefix);    % Magnetometer data

% Load data from MAT files
fprintf('\nLoading sbe file: %s\n', sbe_file);
fprintf('Loading mag file: %s\n', mag_file);

data_sbe = load(sbe_file, 'renav_sbe49');
data_mag = load(mag_file, 'renav_mag');

% Initialize structured data for easier downstream usage
[magData, sbeData] = Mag_initialize_data(data_mag.renav_mag, data_sbe.renav_sbe49);

clear data_sbe data_mag sbe_file mag_file

%%
% -------------------------------------------------------------------------
% Step 2: Spin fitting

fprintf('\n================ Spin fitting calibration =================\n');

spinFitFile  = sprintf('spin_fit_m_%s.mat', file_prefix);   % calibration matrix
spinTimeFile = sprintf('spin_time_m_%s.mat', file_prefix);  % start/end times of spins

if isfile(spinFitFile)
    fprintf('Spin fit for this dive exists. Recalculate? (y/n): ');
    answ = input('', 's');
    if strcmpi(answ, 'y')
        fprintf('\nRecomputing spin fit coefficients...\n');
        Fitting = Mag_spin_fitting(magData, sbeData, spinFitFile, spinTimeFile, sentry_number);
    else
        fprintf('Loading existing spin fit coefficients from file.\n');
        load(spinFitFile, 'Fitting');
    end
else
    fprintf('Computing spin fit coefficients (first time)...\n');
    Fitting = Mag_spin_fitting(magData, sbeData, spinFitFile, spinTimeFile, sentry_number);
end

clear answ

%%
% -------------------------------------------------------------------------
% Step 3: Data processing options

fprintf('\n=================== Processing options ====================\n');
fprintf('\n');
fprintf('1: Process a short-duration segment for a quick validation.\n')
fprintf('2: Process & save entire dive including spins (.mat/.txt).\n')
fprintf('3: Process & save dive excluding spins at start & end.\n')

fprintf('\n');
choice = input('Process? (1/2/3): ', 's');

% Output file names
output_full = sprintf('%s_dive_withspins.mat', file_prefix);
output_clean = sprintf('%s_dive_nospin.mat', file_prefix);

% Process based on user choice
if strcmpi(choice, '1')
    repeat_seg = true;
    while repeat_seg
        fprintf('\nChoose a segment of different duration...\n');
        Mag_process_segment(magData, sbeData, Fitting);
        fprintf('\n#-------------- Segment verification is done! ------------#\n');

        repeat_seg = Mag_ask_yesno('\nRun another segment check? (y/n): ');
    end

    fprintf('\nSegment check complete.\n');

    % Automatically prompt for full processing (option 2 or 3)
    fprintf('\nNow process and save the full dive?\n');
    fprintf('2: Process & save entire dive including spins (.mat/.txt)\n');
    fprintf('3: Process & save dive excluding spins at start & end\n');

    next_choice = input('Choose 2/3: ', 's');

    if strcmpi(next_choice, '2')
        Mag_save_withspins(magData, sbeData, Fitting, file_prefix, output_full);
        
        fprintf('\n---------------------------------------------------------- \n')
        if Mag_ask_yesno('Also save part of dive excluding spins time? (y/n): ')
            Mag_save_nospins(magData, sbeData, Fitting, file_prefix, output_clean);
        end

    elseif strcmpi(next_choice, '3')
        Mag_save_nospins(magData, sbeData, Fitting, file_prefix, output_clean);

        fprintf('\n---------------------------------------------------------- \n')
        if Mag_ask_yesno('Also want to process full dive with spins? (y/n): ')
            Mag_save_withspins(magData, sbeData, Fitting, file_prefix, output_full);
        end
    else
        fprintf('[WARN] Invalid follow-up choice. No processing performed.\n');
    end

elseif strcmpi(choice, '2')
    Mag_save_withspins(magData, sbeData, Fitting, file_prefix, output_full);

    fprintf('\n-------------------------------------------------------------------- \n')
    if Mag_ask_yesno('Also save part of dive excluding spins time? (y/n): ')
        Mag_save_nospins(magData, sbeData, Fitting, file_prefix, output_clean);
    end

elseif strcmpi(choice, '3')
    Mag_save_nospins(magData, sbeData, Fitting, file_prefix, output_clean);

    fprintf('\n-------------------------------------------------------------------- \n')
    if Mag_ask_yesno('Also want to process full dive with spins? (y/n): ')
        Mag_save_withspins(magData, sbeData, Fitting, file_prefix, output_full);
    end

else
    fprintf('[WARN] Invalid choice. No processing performed.\n');
end

clear choice next_choice repeat_seg spinTimeFile spinFitFile Fitting
%%
% -------------------------------------------------------------------------
% Step 4: Crossover correction

fprintf('\n================== Crossover Correction ===================\n');

if Mag_ask_yesno('Perform crossover correction? (y/n): ')
    close all;
    seg_length = input('Enter segment length (e.g., 300): ');
    hdg_thresh = input('Enter heading threshold in degrees (e.g., 50): ');

    % Ask user whether to run only on this dive or multiple dives
    fprintf('\nChoose crossover mode:\n');
    fprintf('1: Perform crossover correction ONLY on this dive (%s)\n', sentry_number);
    fprintf('2: Perform crossover correction ACROSS multiple dives\n');
    mode_choice = input('Select option (1/2): ', 's');

    if strcmp(mode_choice, '2')
        sentryNums = input('Enter sentry number(s) as array (e.g., [744 745]): ');
    else
        sentryNums = str2double(sentry_number);  % convert single input
    end

    for i = 1:length(sentryNums)
        sn = sentryNums(i);

        % Step 4a: Run crossover correction
        fprintf('\nRunning crossover correction for Sentry %d...\n', sn);
        Mag_crossover_correction(seg_length, hdg_thresh, sn);
    end

else
    fprintf('Skipped crossover + anomaly correction.\n');
end

clear folderPath i mode_choice sentryNums sn

fprintf('\n#-------------- Magnetic data processing done! ------------#\n');
function Mag_save_nospins(magData, sbeData, Fitting, file_prefix, output_file)
    % Processes and optionally saves dive data excluding spin times
    % - Checks if the specified output file already exists
    % - If it does, prompts the user whether to recompute the data
    % - If not or if the user chooses to recompute, processes dive data 
    %   while excluding spin times
    %
    % Syntax:
    %   save_nospins(magData, sbeData, Fitting, file_prefix, output_file)
    %
    % Inputs:
    %   magData      : Struct or matrix of magnetometer data used for dive processing
    %   sbeData      : Struct or matrix of SBE (CTD/sensor) data required for dive analysis
    %   Fitting      : Struct of fitting parameters or model coefficients for correction
    %   file_prefix  : String used as a prefix for any intermediate or output files
    %   output_file  : Full path to the output file for storing processed data

%-------------------------------------------------------------------------
% Step 1: Check if file exists or not and get the required data
    if isfile(output_file)
        if Mag_ask_yesno(sprintf('File "%s" exists. Recompute? (y/n): ', output_file))
            fprintf('Recomputing no-spin output...\n');
            Mag_process_dive_without_spins(magData, sbeData, Fitting, file_prefix);
        else
            fprintf('Skipped recomputation of dive without spins.\n');
        end
    else
        fprintf('\nProcessing and saving dive exluding spins time...\n');
        Mag_process_dive_without_spins(magData, sbeData, Fitting, file_prefix);
    end
end
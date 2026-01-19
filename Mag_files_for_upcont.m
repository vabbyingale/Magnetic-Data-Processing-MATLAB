clear all; close all; clc

% === File selection ===
allFiles = dir('Seg_sentry582*_cross.txt');
pattern = '^Seg_sentry\d{3}_\d{8}_\d{4}_cross\.txt$';

fileList = {};
for k = 1:length(allFiles)
    fileName = allFiles(k).name;
    if ~isempty(regexp(fileName, pattern, 'once'))
        fileList{end+1} = fileName;
    end
end
fileList = sort(fileList);
numDives = length(fileList);

% === Loop over files ===
for i = 1:numDives
    inFile = fileList{i};
    outFile = strrep(inFile, '_cross.txt', '_corrected.txt');
    
    % Load file
    data = readmatrix(inFile);

    % Remove rows with any NaN values
    data = data(~any(isnan(data),2),:);

    % Coordinates
    lat = data(:,2); 
    lon = data(:,3);

    % Use first valid index
    lat1 = lat(find(~isnan(lat), 1));  
    lon1 = lon(find(~isnan(lon), 1));

    % Declination and Inclination (deg)
    Dec = 11.0707; 
    Inc = -21.5330;

    % IGRF values (N, E, D)
    J = Mag_magfd(2024.9, 1, 0.003, 90 - lat1, lon1);  
    IGRF = [J(1), J(2), J(3)];

    % Anomalies (subtract IGRF)
    anom_n = data(:,25) - IGRF(1);   % North anomaly
    anom_e = data(:,26) - IGRF(2);   % East anomaly
    anom_d = data(:,27) - IGRF(3);   % Down anomaly

    % Total-field anomaly (projected along field direction)
    anom_t = anom_n .* cosd(Dec) .* cosd(Inc) + ...
             anom_e .* sind(Dec) .* cosd(Inc) + ...
             anom_d .* sind(Inc);

    % Append anomalies as new columns
    data_corrected = [data anom_n anom_e anom_d anom_t];

    % Save corrected file (comma delimiter, 6 decimals)
    dlmwrite(outFile, data_corrected, 'delimiter', ',', 'precision', '%.6f');

end

disp(['Processed ' num2str(numDives) ' files successfully.']);

function Mag_crossover_correction(segment_length, heading_threshold, sentryNums)
% MAG_CROSSOVER_CORRECTION_ANOM
% Perform crossover correction on Sentry magnetic anomaly data
%
% T = Mag_crossover_correction_anom(segment_length, heading_threshold, sentryNums)
%
% Inputs:
%   segment_length    - Segment length (e.g., 300)
%   heading_threshold - Heading difference threshold (degrees, e.g., 50)
%   sentryNums        - One or more sentry numbers (e.g., [744 745])
%
% Output:
%   T - Table with RMS before and after correction for each magnetic variable
%

    if nargin < 3 || isempty(sentryNums)
        sentryNums = input('Enter sentry number(s) (e.g., 582 or [582 583]): ');
    end

    % === Discover all dive files ===
    allFiles = [];
    for s = sentryNums
        patternStr = sprintf('Seg_sentry%03d*.txt', s);
        allFiles = [allFiles; dir(patternStr)]; %#ok<AGROW>
    end

    pattern = '^Seg_sentry\d{3}_\d{8}_\d{4}\.txt$';
    fileList = {};
    for k = 1:length(allFiles)
        if ~isempty(regexp(allFiles(k).name, pattern, 'once'))
            fileList{end+1} = allFiles(k).name; %#ok<AGROW>
        end
    end
    fileList = sort(fileList);

    numDives = length(fileList);
    if numDives == 0
        error('No valid sentry files found.');
    end

    % === Magnetic variables ===
    magVars = {'Anom_n','Anom_e','Anom_d','Anom_total'};
    rmsBeforeAll = zeros(numel(magVars),1);
    rmsAfterAll  = zeros(numel(magVars),1);

    % === Load all dive files once ===
    diveData = struct;
    for d = 1:numDives
        data = readtable(fileList{d});
        diveData(d).raw = data;
        diveData(d).numSegments = floor(height(data)/segment_length);
    end

    % === Loop over magnetic variables ===
    for varIdx = 1:numel(magVars)
        magVar = magVars{varIdx};
        fprintf('Processing %s: ', magVar);

        % === Segment all dives ===
        segment = struct;
        segIdx = 1;

        for d = 1:numDives
            data = diveData(d).raw;

            lat = data.lat;
            lon = data.lon;
            mag = data.(magVar);
            hdg = data.Heading;

            N = diveData(d).numSegments;

            for i = 1:N
                idx = (1:segment_length) + (i-1)*segment_length;

                segment(segIdx).lat      = lat(idx);
                segment(segIdx).lon      = lon(idx);
                segment(segIdx).mag      = mag(idx);
                segment(segIdx).hdg      = hdg(idx);
                segment(segIdx).dive     = d;
                segment(segIdx).localIdx = i;

                segIdx = segIdx + 1;
            end
        end

        numSegments = numel(segment);

        % === Find crossovers ===
        crossoverData = [];

        for i = 1:numSegments-1
            for j = i+1:numSegments

                if (segment(i).dive == segment(j).dive) && ...
                   (abs(segment(i).localIdx - segment(j).localIdx) <= 1)
                    continue
                end

                for m = 1:numel(segment(i).lat)

                    latA = segment(i).lat(m);
                    lonA = segment(i).lon(m);
                    magA = segment(i).mag(m);
                    hdgA = segment(i).hdg(m);

                    dists = haversine(latA, lonA, segment(j).lat, segment(j).lon);
                    [minDist, idxMin] = min(dists);

                    if minDist < 0.01  % km
                        magB = segment(j).mag(idxMin);
                        hdgB = segment(j).hdg(idxMin);

                        dh = abs(hdgA - hdgB);
                        dh = min(mod(dh,360), 360 - mod(dh,360));

                        if dh > heading_threshold
                            diffMag = magA - magB;
                            crossoverData = [crossoverData; i, j, diffMag]; %#ok<AGROW>
                        end
                    end
                end
            end
        end

        % === Least Squares Correction ===
        A = zeros(size(crossoverData,1), numSegments);
        b = zeros(size(crossoverData,1),1);

        for k = 1:size(crossoverData,1)
            A(k, crossoverData(k,1)) =  1;
            A(k, crossoverData(k,2)) = -1;
            b(k) = -crossoverData(k,3);
        end

        x = A \ b;

        % === Apply corrections to segments ===
        for i = 1:numSegments
            segment(i).mag_corrected = segment(i).mag + x(i);
        end

        % === Reconstruct corrected dives ===
        for d = 1:numDives
            data = diveData(d).raw;
            N = diveData(d).numSegments;

            mag_corrected = NaN(height(data),1);

            for i = 1:N
                idx = (1:segment_length) + (i-1)*segment_length;
                gIdx = find([segment.dive]==d & [segment.localIdx]==i);
                mag_corrected(idx) = segment(gIdx).mag_corrected;
            end

            correctedVarName = sprintf('%s_crossover', magVar);
            diveData(d).raw.(correctedVarName) = mag_corrected;
        end

        % === RMS statistics ===
        diff0 = crossoverData(:,3);
        diff1 = diff0 + x(crossoverData(:,1)) - x(crossoverData(:,2));

        rmsBeforeAll(varIdx) = rms(diff0);
        rmsAfterAll(varIdx)  = rms(diff1);

        fprintf('RMS Before: %.4f nT | RMS After: %.4f nT\n', ...
                rmsBeforeAll(varIdx), rmsAfterAll(varIdx));
    end

    % === Save corrected files ===
    for d = 1:numDives
        data = diveData(d).raw;
        [~, name, ext] = fileparts(fileList{d});
        outFile = [name '_cross' ext];
    
        % Open file for writing
        fid = fopen(outFile, 'w');
    
        % Write headers
        headers = data.Properties.VariableNames;
        fprintf(fid, '%s', headers{1});
        for c = 2:numel(headers)
            fprintf(fid, ',%s', headers{c});
        end
        fprintf(fid, '\n');
    
        % Number of columns
        nCols = width(data);
    
        % Indices of last 4 columns
        last4Idx = (nCols-3):nCols;
    
        % Write data row by row
        for row = 1:height(data)
            for col = 1:nCols
                val = data{row, col};
                if ismissing(val)
                    % Handle missing data (NaN or <missing>)
                    str = '';
                elseif ismember(col, last4Idx)
                    % Last 4 columns with 6 decimals
                    str = sprintf('%.6f', val);
                elseif isnumeric(val)
                    % Other numeric columns default format
                    str = num2str(val);
                elseif isstring(val) || ischar(val) || iscellstr(val)
                    % String columns (if any)
                    str = char(val);
                else
                    % For other data types fallback to string conversion
                    str = string(val);
                end
    
                if col == 1
                    fprintf(fid, '%s', str);
                else
                    fprintf(fid, ',%s', str);
                end
            end
            fprintf(fid, '\n');
        end
    
        fclose(fid);
    end

    fprintf('Saved crossover corrected file: %s\n', outFile);

end

% === Haversine distance function ===
function d = haversine(lat1, lon1, lat2, lon2)
    R = 6371; % km
    dlat = deg2rad(lat2 - lat1);
    dlon = deg2rad(lon2 - lon1);
    a = sin(dlat/2).^2 + ...
        cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dlon/2).^2;
    d = 2 * R * atan2(sqrt(a), sqrt(1-a));
end

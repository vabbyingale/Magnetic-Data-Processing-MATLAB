function [seg_start, seg_end] = Mag_select_dive_segment(sbe_heading, sbe_lon, sbe_lat, sbe_utime, sbe_alt)
    % A function identify and select a dive track that also includes a spin maneuver.
    %
    % Syntax:
    %   [start, end] = dive_segment_entire(heading, lon, lat, utime, alt);
    %
    % Inputs:
    %   sbe_heading - Vector of heading angles (degrees)
    %   sbe_lon     - Vector of longitudes
    %   sbe_lat     - Vector of latitudes
    %   sbe_utime   - Vector of Unix times
    %   sbe_alt     - Vector of altitudes (m)
    %
    % Outputs:
    %   seg_start   - Unix time marking start of selected track
    %   seg_end     - Unix time marking end of selected track


    %-------------------------------------------------------------------------
    % Step 1: Identify heading jumps

    ang = diff(sbe_heading);             % get heading difference
    ang = [ang;0];                       % pad to match original length
    ndc = find(abs(ang)>180);            % indices where heading exceeds > 180°
    jumps = ang(ndc);                    % store jump values for reference

    %-------------------------------------------------------------------------
    % Step 2: Plotting

    % Plot track with spin indicators
    figure(1);        
    plot(sbe_lon, sbe_lat, 'b-', 'LineWidth', 0.8, 'DisplayName', 'Track'); hold on
    plot(sbe_lon(ndc), sbe_lat(ndc), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Heading change'); legend
    hold off;

    % Plot altitude vs time with heading indicators
    figure(2);                 
    plot(sbe_utime, sbe_alt, 'b-', 'LineWidth', 0.8); hold on
    plot(sbe_utime(ndc), sbe_alt(ndc), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6)
    hold off;

    % Plot heading vs time with heading indicators
    figure(3);                 
    plot(sbe_utime, sbe_heading, 'b-', 'LineWidth', 0.8); hold on
    plot(sbe_utime(ndc), sbe_heading(ndc), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6)
    hold off;

    %-------------------------------------------------------------------------
    % Step 3: First zoom and first selection

    format longG;    
    fprintf('\n');
    disp('Zoom in Figure 3 for first limit')
    disp('When ready, press any key to confirm zoom');
    zoom on;
    pause;      % wait for any key press to confirm zoom
    zoom off;

    figure(3);
    fprintf('\n')
    disp('Now click to select the first limit on the plot.');
    [x1, ~] = ginput(1);

    %-------------------------------------------------------------------------
    % Step 4: Second zoom and second selection

    disp('Change the zoom to adjust second limit') 
    disp('When ready, press any key to confirm zoom.');
    zoom on;
    pause;
    zoom off;

    figure(3);
    fprintf('\n')
    disp('Now click to select the second limit on the plot.');
    [x2, ~] = ginput(1);

    %-------------------------------------------------------------------------
    % Step 5: Arrange limits and zoom plots accordingly

    seg_start = min(x1, x2);
    seg_end   = max(x1, x2);

    figure(3);
    xlim([seg_start seg_end]);

    figure(2);
    xlim([seg_start seg_end]);

    %-------------------------------------------------------------------------
    % Step 6: Plot dive track segment in new figure

    figure(4);
    track_ndc = find(sbe_utime > seg_start & sbe_utime < seg_end);
    plot(sbe_lon, sbe_lat, 'b-', 'LineWidth', 0.8); hold on
    plot(sbe_lon(track_ndc), sbe_lat(track_ndc), 'r-', 'LineWidth', 1.5);
    xlabel('Longitude'); ylabel('Latitude');
    title('Selected Dive Track Segment');
    hold off;
end

function [spin_start, spin_end] = Mag_find_spins(sbe_heading, sbe_lon, sbe_lat, sbe_utime, sbe_alt, sentry_number)
    % A function to identify and select a time segment corresponding to only spin maneuver
    %
    % Syntax:
    %   [start, end] = find_spins(heading, lon, lat, utime, alt);
    %
    % Inputs:
    %   sbe_heading = Vector of heading angles (degrees)
    %   sbe_lon     = Vector of longitudes
    %   sbe_lat     = Vector of latitudes
    %   sbe_utime   = Vector of Unix times
    %   sbe_alt     = Vector of altitudes (meters)
    %
    % Outputs:
    %   spin_start  = Unix time marking start of selected spin
    %   spin_end    = Unix time marking end of selected spin


%-------------------------------------------------------------------------
% Step 1: Identify heading jumps
ang = diff(sbe_heading);             % get heading difference
ang = [ang;0];                       % pad to match original length
ndc = find(abs(ang)>180);            % indices where heading exceeds > 180°
jumps = ang(ndc);                    % store jump values for reference

%-------------------------------------------------------------------------
% Step 2: Plotting

% Plot track with spin indicators
fig1 = figure('Position',[200, 200, 600, 500]);        
plot(sbe_lon, sbe_lat, 'k-', 'LineWidth', 0.8,'DisplayName','Track'); hold on
plot(sbe_lon(ndc), sbe_lat(ndc), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Heading change'); 
legend('Location','southeast')
xlabel('Longitude'); ylabel('Latitude'); grid on
Mag_position_plot(gca)
Mag_save_plot_prompt(fig1, sprintf('01_S%s_track', sentry_number));

% Plot altitude vs time with heading indicators
figure(2)                 
plot(sbe_utime, sbe_alt, 'k-', 'LineWidth', 0.8); hold on; grid on
plot(sbe_utime(ndc), sbe_alt(ndc), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize',6)

% Plot heading vs time with heading indicators
figure(3)                  
plot(sbe_utime, sbe_heading, 'k-', 'LineWidth', 0.8); hold on; grid on
plot(sbe_utime(ndc), sbe_heading(ndc), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6)

%-------------------------------------------------------------------------
% Step 3: Ask user to select start and end times
format longG;   

disp('Zoom in Figure 3 near spins....')
disp('Press any key & select the spin limits by click!') 
fprintf('\n')
pause;

% Select two points on heading plot
[x,y] = ginput(1);      p1 = [x(1),y(1)];
[x,y] = ginput(1);      p2 = [x(1),y(1)];

% Ensure seg_start < seg_end
if p1(1) > p2(1)                  
    spin_start = p2(1);
    spin_end = p1(1);
else                              
    spin_start = p1(1);
    spin_end = p2(1);
end

% Zoom in to selected time window on heading plot
f1 = figure(3);
xlim([spin_start spin_end])
ylabel('Heading (°)')
xlabel('Unix Time')
Mag_position_plot(gca)
Mag_save_plot_prompt(f1, sprintf('02_S%s_spin_heading', sentry_number));

% Zoom in to selected time window on altitude plot
f2 = figure(2);
xlim([spin_start spin_end])
ylabel('Altitude (m)')
xlabel('Unix Time')
Mag_position_plot(gca)
Mag_save_plot_prompt(f2, sprintf('03_S%s_spin_altitude', sentry_number));

% Plot spins in new figure
f3 = figure(4);
track_ndc = find(sbe_utime > spin_start & sbe_utime < spin_end);
plot(sbe_lon(track_ndc), sbe_lat(track_ndc), 'k-', 'LineWidth', 0.8); grid on
Mag_position_plot(gca)
Mag_save_plot_prompt(f3, sprintf('04_S%s_spins', sentry_number));

hold off;
end



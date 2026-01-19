function [start_utime, end_utime] = Mag_find_straight_segment(duration,sbe_utime,sbe_lon,sbe_lat,mag_utime)
    % A function that extracts straight segment of dive track
    %
    % Syntax:
    %   [start_utime, end_utime] = find_straight_segment(120, sbe_utime, sbe_lon, sbe_lat, mag_utime);
    %
    % Inputs:
    %   duration   = Duration of straight segment to extract (in seconds)
    %   sbe_utime  = Unix time vector
    %   sbe_lon    = Longitude from Sentry navigation
    %   sbe_lat    = Latitude from Sentry navigation
    %   mag_utime  = Unix time vector from Maggie (magnetometer) file, used to determine sampling rate
    %
    % Outputs:
    %   start_utime = Start Unix time of selected straight segment
    %   end_utime   = End Unix time of selected straight segment

%-------------------------------------------------------------------------
% Step 1: Compute sampling rate of magnetometer data
mag_samp_rate = round(1/((mag_utime(end)-mag_utime(1))/length(mag_utime)));     % Sampling rate 
mag_nsamp = round(duration/mag_samp_rate);                                      % Number of samples 
sbe_utime_mid = (sbe_utime(length(sbe_utime))-sbe_utime(1))/2 + sbe_utime(1);   % find midpoint time of dive

%-------------------------------------------------------------------------
% Step 2: calculate R2 (linear regression)
rsq_max = 0;                                                        % Max R2 value
i_max = 0;                                                          % Index for max R2 value

for i = 1:3                                    
 FG_ind = find(abs(mag_utime-sbe_utime_mid)<duration/2);
 rsq = fitlm(sbe_lon(FG_ind),sbe_lat(FG_ind)).Rsquared.Ordinary;    % Compute R2

 % Keep segment with highest R2
 if rsq > rsq_max                        
     i_max = i;
     rsq_max = rsq;
     sbe_utime_best = sbe_utime_mid;                              
 end
 sbe_utime_mid = sbe_utime_mid + duration;                          % move window forward for iteration
end

%-------------------------------------------------------------------------
% Step 3: Plotting a segment
figure('Position',[100, 100, 600, 600]);
plot(sbe_lon, sbe_lat, 'b-', 'LineWidth', 0.8); hold on
mag_ind = find(abs(sbe_utime - sbe_utime_best) < duration/2);       % Get best indices from mag file
plot(sbe_lon(mag_ind), sbe_lat(mag_ind), 'r-', 'LineWidth', 3)      % Plot straight line segment
xlabel("Longitude"); ylabel("Latitude")
title(sprintf('Best Straight Segment (%.0f s) - R^2 = %.3f', duration, rsq_max));
legend({'Full Dive Track', 'Straight Segment'}, 'Location', 'best');

ax1 = gca; 
Mag_position_plot(ax1)
hold off;

% -------------------------------------------------
% Step 4: compute start and end Unix times for segment
start_utime = sbe_utime_best-duration/2;
end_utime = sbe_utime_best+duration/2;

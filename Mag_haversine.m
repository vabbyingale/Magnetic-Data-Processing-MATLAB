function d = Mag_haversine(lat1, lon1, lat2, lon2)
    % Compute great-circle distance between two coordinates
    %
    % Syntax:
    %   d = haversine(lat1, lon1, lat2, lon2)
    %
    % Inputs:
    %   lat1, lon1 - Latitude and longitude of the first point (degrees)
    %   lat2, lon2 - Latitude and longitude of the second point (degrees)
    %
    % Output:
    %   d - Distance between the two points along the surface of the Earth (km)

% -------------------------------------------------------------------------
% Step 1: Convert input coordinates from degrees to radians
R = 6371;                                     % Earth's radius in km
dlat = deg2rad(lat2 - lat1);                  % change in lat
dlon = deg2rad(lon2 - lon1);                  % change in lon

% -------------------------------------------------------------------------
% Step 2: Apply the haversine formula
a = sin(dlat/2).^2 + cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));
d = R * c;                                    % Distance in km

function [mag_ned] = Mag_component_rotation(heading, pitch, roll, XYZ_tcal, orientation)
    % A function to rotate calibrated XYZ magnetic field measurements into NED frame.  
    %
    % Syntax:
    %   mag_ned = ComponentRotation(heading, pitch, roll, XYZ_tcal, orientation)
    %
    % Inputs:
    %   heading     = Vector of heading (degrees)
    %   pitch       = Vector of pitch (degrees)
    %   roll        = Vector of roll (degrees)
    %   XYZ_tcal    = Nx3 matrix of temperature-corrected magnetic field components (X,Y,Z)
    %   orientation = Integer flag specifying sensor mounting orientation
    %                 0: X =  X, Y =  Y, Z =  Z
    %                 1: X = -X, Y = -Y, Z =  Z (top mount)
    %                 2: X = -X, Y =  Y, Z = -Z (port mount)
    %                 3: X =  X, Y = -Y, Z = -Z (starboard mount)
    %                 4: X = -X, Y =  Y, Z =  Z
    %                 5: X =  X, Y =  Y, Z = -Z
    %
    % Outputs:
    %   mag_ned     = Nx3 matrix of rotated magnetic field components in NED frame (North, East, Down)
 
% -------------------------------------------------------------------------
% Step 1: prepare input matrix and output vector

% Convert degrees to radians
alpha=heading*pi/180;               % Heading
beta=pitch*pi/180;                  % Pitch
gamma=roll*pi/180;                  % Roll

mag_ned=zeros(size(XYZ_tcal));      % output ned components

% -------------------------------------------------------------------------
% Step 2: Adjust sensor axes based on mounting orientation

% Axis orientation for different configuration
switch orientation

    case 0                        % No changes needed
        XYZ_tcal(:,1)=XYZ_tcal(:,1);      % X
        XYZ_tcal(:,2)=XYZ_tcal(:,2);      % Y
        XYZ_tcal(:,3)=XYZ_tcal(:,3);      % Z

    case 1                        % Top mount or FG2: (X=-X, Y=-Y, Z=Z)
        fprintf('\n')
        disp('Component rotation for top maggie')
        XYZ_tcal(:,1)=-XYZ_tcal(:,1);      % X
        XYZ_tcal(:,2)=-XYZ_tcal(:,2);      % Y

    case 2                        % Port mount or FG1: (X=-X, Y=Y, Z=-Z)
        fprintf('\n')
        disp('Component rotation for port maggie')
        XYZ_tcal(:,1) = -XYZ_tcal(:,1);    % X
        XYZ_tcal(:,3) = -XYZ_tcal(:,3);    % Z

    case 3                        % Starboard or FG0: (X=X, Y=-Y, Z=-Z)
        fprintf('\n')
        disp('Component rotation for starboard maggie')
        XYZ_tcal(:,2) = -XYZ_tcal(:,2);
        XYZ_tcal(:,3) = -XYZ_tcal(:,3); 

    case 4                        % X = -X
        XYZ_tcal(:,1) = -XYZ_tcal(:,1);

    case 5                        % Z = -Z
        XYZ_tcal(:,3) = -XYZ_tcal(:,3);

    otherwise
        error('Unknown orientation code: %d', orientation);
end

% -------------------------------------------------------------------------
% Step 3: Rotate each sample to NED

% precompute trigemometric function
for k = 1:min(length(heading), size(XYZ_tcal, 1))
    sa=sin(alpha(k));   sb=sin(beta(k));    sg=sin(gamma(k));
    ca=cos(alpha(k));   cb=cos(beta(k));    cg=cos(gamma(k));

% build rotation matrix T
    T(1,1)=ca*cb;       T(2,1)=ca*sb*sg-sa*cg;      T(3,1)=ca*sb*cg+sa*sg;
    T(1,2)=sa*cb;       T(2,2)=sa*sb*sg+ca*cg;      T(3,2)=sa*sb*cg-ca*sg;
    T(1,3)=-sb;         T(2,3)=cb*sg;               T(3,3)=cb*cg;
    
% Rotate the calibrated XYZ vector into NED
    mag_ned(k,:) = (T'*XYZ_tcal(k,:)')';
end
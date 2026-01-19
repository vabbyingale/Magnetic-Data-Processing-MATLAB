function [COR,fcorT] = Mag_calib1520_sensors(X, sensor)
    % A function that applies temperature and orientation calibration to 3-axis magnetometer data
    %
    % Syntax:
    %   [COR, fcorT] = calib1520_sensors(X, sensor)
    %
    % Inputs:
    %      X   = Nx4 matrix of raw measurements [X Y Z Temp], in gauss and °C
    % sensor   = Sensor serial number (e.g. 0687, 0688, 0689, 0690) to select calibration curves
    %
    % Outputs:
    %   COR    = Nx3 corrected field components (nT)
    %   fcorT  = 1xN total corrected field strength (nT)

% -------------------------------------------------------------------------
% Step 1: Convert input gauss measurements to nT (1 gauss = 1e5 nT)
raw(:,1)=X(:,1) * 1e5;          % X component
raw(:,2)=X(:,2) * 1e5;          % Y component
raw(:,3)=X(:,3) * 1e5;          % Z component
T=X(:,4);                       % Temperature

fprintf('\n')
fprintf('For temperature correction:\n')
disp('choose sensor 688: starboard / 689: port / 690: top')
sensor = input('Enter sensor number: ');

% -------------------------------------------------------------------------
% Step 2: Select sensor-specific calibration coefficients
% Polynomial bias corrections (offsets), scale factors, rotations
if sensor == 0687
    b(:,1) = 25.837 + 1.3642 .* T - 0.058418 .* T.^2;                 
    b(:,2) = 84.498 - 4.6222 .* T + 0.048285 .* T.^2;
    b(:,3) = 37.862 - 2.8216 .* T + 0.0032608 .* T.^2;

    scale(:,1) = 0.99397 + 2.1563e-4 .* T + 4.1443e-7 * T.^2;
    scale(:,2) = 0.99141 + 2.6452e-4 .* T + 5.8845e-6 * T.^2;
    scale(:,3) = 0.99623 + 1.7733e-4 .* T + 1.6337e-6 * T.^2 ;

    r(:,1) = -197.81 + 33.499 .* T - 0.69658 * T.^2 ;
    r(:,2) = 375.91 + 17.633 .* T - 0.46883 * T.^2 ;
    r(:,3) = 238.68 - 11.479 .* T + 0.044622 * T.^2 ;

elseif sensor == 0688
    b(:,1) = 108.28 - 4.4685 .* T + 1.5674e-3 .* T.^2;
    b(:,2) = 31.393 - 0.72391 .* T - 3.0317e-2 .* T.^2;
    b(:,3) = -45.188 + 1.6138 .* T - 9.6169e-3 .* T.^2;

    scale(:,1) = 0.99784 + 6.976e-6 .* T + 6.4954e-6 * T.^2;
    scale(:,2) = 0.98097 + 4.6894e-4 .* T - 5.682e-6 * T.^2;
    scale(:,3) = 0.99800 + 2.9658e-5 .* T + 5.8352e-6 * T.^2;

    r(:,1) = -199.24 + 11.182 .* T - 0.42256 .* T.^2;
    r(:,2) = 282.09 - 8.9838 .* T - 0.084534 .* T.^2;
    r(:,3) = 42.962 - 20.321 .* T+ 0.28916 .* T.^2;

elseif sensor == 0689
    b(:,1) = 32.603 - 1.7135 .* T + 6.7599e-2 .* T.^2;
    b(:,2) = -62.389 + 7.4899 .* T - 1.5924e-1 .* T.^2;
    b(:,3) = 42.883 - 1.2902 .* T - 5.5892e-2 .* T.^2;
    
    scale(:,1) = 0.99537 + 1.142e-4 .* T + 7.7696e-6 .* T.^2;
    scale(:,2) = 0.99078 + 2.468e-4 .* T + 7.1949e-6 .* T.^2 ;
    scale(:,3) = 0.99659 + 1.331e-5 .* T + 8.2567e-6 .* T.^2;

    r(:,1) = 113.65 + 1.7461 .* T - 0.28297 .* T.^2;
    r(:,2) = -475.47 + 6.822 .* T - 0.41443 .* T.^2;
    r(:,3) = 91.459 - 12.848 .* T + 0.53951 .* T.^2;

elseif sensor == 0690
    b(:,1) = -58.421 - 1.958 .* T + 6.1023e-2 .* T.^2;
    b(:,2) = 41.01 - 8.5493 .* T + 3.0701e-1 .* T.^2;
    b(:,3) = 227.85 + 4.0855 .* T - 1.8867e-1 .* T.^2 ; 

    scale(:,1) = 0.99381 + 2.4977e-4 .* T + 2.7546e-6 .* T.^2;
    scale(:,2) = 0.98771 + 6.3261e-4 .* T - 5.0385e-6 .* T.^2;
    scale(:,3) = 0.99634 + 1.2419e-4 .* T + 5.662e-6 .* T.^2;

    r(:,1) = 435.02 - 9.6061 .* T+ 0.20449 .* T.^2;
    r(:,2) = 193.36 + 6.2856 .* T + 0.1021 .* T.^2;
    r(:,3) = -189.36 + 2.2936 .* T - 0.18431 .* T.^2;        
end

% -------------------------------------------------------------------------
% Step 3: Initialize outputs
fraw = zeros(length(raw),1);        % uncorrected field magnitude
fcor = zeros(length(raw),1);        % corrected field magnitude
COR  = zeros(length(raw),3);        % corrected XYZ matrix

% -------------------------------------------------------------------------
% Step 4: Apply scale, bias, and orientation correction for each sample
for k=1:length(raw)
S=diag(scale(k,:));                 % scale matrix (diagonal)
Sinv=inv(S);                        % inverse matrix

u=(r(k,:)/3600)*pi/180;             % rotation angles (r) in arcseconds, convert to radians
su=sin(u); cu=cos(u);               % trig components for rotation matrix
w=sqrt(1-su(2)^2-su(3)^2);          % ensure rotation matrix is normalized

% construct 3x3 transformation matrix P
P=zeros(3);
P(1,1)=1;   P(2,1)=su(1)/cu(1);     P(2,2)=1/cu(1);
P(3,1)=-((su(1)*su(3)+cu(1)*su(2))/(w*cu(1)));
P(3,2)=-su(3)/(w*cu(1));
P(3,3)=1/w;

PS=P*Sinv;                          % combine P and invserse scale matrix
tmp=raw(k,:)-b(k,:);                % apply bias correction

% apply full correction (scale + rotation) 
cor=zeros(1,3);
for i=1:3
    for j=1:3
        cor(i)=cor(i)+PS(i,j)*tmp(j);
    end
end

% Store uncorrected and corrected field magnitudes
fraw(k)=sqrt(raw(k,1)^2+raw(k,2)^2+raw(k,3)^2);
fcor(k)=sqrt(cor(1)^2+cor(2)^2+cor(3)^2);
COR(k,:)=cor;
end

% -------------------------------------------------------------------------
% Step 5: Return corrected field vectors and transposed magnitude
fcorT=fcor';

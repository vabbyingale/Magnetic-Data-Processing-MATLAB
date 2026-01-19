function  J = Mag_magfd(DATE,ITYPE,ALT,COLAT,ELONG)
    % MAGFD  Compute Earth's magnetic field using IGRF model using sh2020.mat file 
    % ONLY uses harmonic expansion to order 10 at present.
    %
    % Usage:
    %   J = magfd(DATE, ITYPE, ALT, COLAT, ELONG)
    %
    % Inputs:
    %   DATE  = date of survey (decimal years)
    %   ITYPE = 1 for geodetic coordinates, 2 for geocentric coordinates
    %   ALT   = altitude above sea level in km for ITYPE = 1 and
    %         = radial distance from earth's center if ITYPE = 2
    %   COLAT = colatitute (90 - latitude) in decimal degrees
    %   ELONG = longitude in decimal degrees
    %
    % Output:
    %   J     = 1x4 vector in nT: 
    %         = X as north, Y as east, Z as down and T as total
    %
    % Reference:
    %   IAGA, Division V, Working Group VMOD, 
    %   The 10th generation International Geomagnetic 
    %   Reference Field, Geophys. J. Int, 161, 561-565, 2005.
    %
    % Maurice A. Tivey March 1997
    % http://deeptow.whoi.edu/matlab.html
    % Copyright: Maurice A. Tivey, 2005
    % Woods Hole Oceanographic Institution

%-------------------------------------------------------------------------
% Step 1: Demo if no args given
if nargin < 1
 fprintf('=======DEMO========\n')
 disp('Magnetic field at given lat lon using IGRF model')
 disp('Example: field at Woods Hole from Jan 1st 1900');
 disp('Latitude 42 N, Longitude 74W')
 disp('sample command:  out=magfd(1997,1,0,90-42,-74);')

 % Field at WHOI at 5 yr intervals
 for i=1:25
    out(i,:) = magfd(-((i-1)*5+1900),1,0,90-42,-74);
 end
 plot(([1:25]-1)*5+1900, out(:,4), '-r+', 'linewidth', 2);
 xlabel('Year');
 ylabel('Total Magnetic Field (nT)')
 title('Total Magnetic Field Intensity at Woods Hole, MA')
 axis tight
 return
end

%-------------------------------------------------------------------------
% Step 2: Determine base year & load spherical harmonic data

    DGRF = 1000:5:2015;
igrfyear = 2020;
igrffile ='sh2020';
      pl = 0;    
if DATE  < 0, pl = 1; end
    DATE = abs(DATE);

% Determine year for base DGRF to use.
 if DATE < igrfyear
     BASE = fix(DATE-DGRF(1)); 
        i = fix(BASE/5) + 1;
     BASE = DGRF(i);
     
     if pl==0
         fprintf('Using DGRF base year %f \n',BASE);      
     end
     eval(['load sh',num2str(BASE)])
  
     iagh = agh; iagh41 = agh41;       % loads agh and agh41

     if BASE < 1900                    % a check to get pre-1900 estimates of gauss coeffs
        eval(['load sh',num2str(BASE+25)])

     else
        eval(['load sh',num2str(DGRF(i+1))])
     end

     eagh  = agh; eagh41 = agh41;
     dgh   = (eagh-iagh)./5;
     dgh41 = (eagh41-iagh41)./5;
     agh = iagh;agh41=iagh41;
     clear iagh iagh41 eagh eagh41
     T = DATE - BASE;

 else
     eval(['load ',igrffile])   % load in igrf data file
     T = DATE - igrfyear;
 end

%-------------------------------------------------------------------------
% Step 3: Combine spectral harmonic coefficients up to degree 10
agh = [agh,agh41];
dgh = [dgh,dgh41];

% Setup coordinates
D2R   = pi/180;
R     = ALT;
SLAT  = cos(COLAT*D2R);
CLAT  = sin(COLAT*D2R);
CL(1) = cos(ELONG*D2R);
SL(1) = sin(ELONG*D2R);

% Setup variables
X  = 0.0;    Y  = 0.0;    Z  = 0.0;
CD = 1.0;    SD = 0.0;    L  = 1;
M  = 1;       N = 0;      RE = 6371.2; % Earth's mean radius

if ITYPE == 1
    % convert geodetic to geocentric coordinates using WGS84
	A2    = 6378.137^2;         % squared semi major axis
	B2    = 6356.7523142^2;     % squared semi minor axis
	ONE   = A2*CLAT*CLAT;
	TWO   = B2*SLAT*SLAT;
	THREE = ONE + TWO;
	FOUR  = sqrt(THREE);
	R     = sqrt(ALT*(ALT + 2.0*FOUR) + (A2*ONE + B2*TWO)/THREE);
	CD    = (ALT + FOUR)/R;
	SD    = (A2 - B2)/FOUR*SLAT*CLAT/R;
	ONE   = SLAT;
	SLAT  = SLAT*CD - CLAT*SD;
	CLAT  = CLAT*CD +  ONE*SD;
end

RATIO = RE/R;                   % if geocentric coordinates desired then only define

%-------------------------------------------------------------------------
% Step 4: compute Schmidt quasi-normalized coefficients
P(1)  = 2.0*SLAT;
P(2)  = 2.0*CLAT;
P(3)  = 4.5*SLAT*SLAT - 1.5;
P(4)  = sqrt(27)*CLAT*SLAT;
Q(1)  = -CLAT;
Q(2)  =  SLAT;
Q(3)  = -3.0*CLAT*SLAT;
Q(4)  = sqrt(3)*(SLAT*SLAT - CLAT*CLAT);

%-------------------------------------------------------------------------
% Step 5: loop over harmonic series
 NMAX = 10;                      % Max number of harmonic degrees
 NPQ  = (NMAX * (NMAX + 3))/2;
for K = 1:NPQ
    if N < M 
        M  = 0;
        N  = N + 1;
        RR = RATIO^(N + 2);
        FN = N;
    end
    FM = M;

    if K  >= 5 %8,5,5
        if (M-N) == 0 %,7,6,7
              ONE = sqrt(1.0 - 0.5/FM);
		        J = K - N - 1;
		     P(K) = (1.0 + 1.0/FM)*ONE*CLAT*P(J);
		     Q(K) = ONE*(CLAT*Q(J) + SLAT/FM*P(J));
		    SL(M) = SL(M-1)*CL(1) + CL(M-1)*SL(1);
		    CL(M) = CL(M-1)*CL(1) - SL(M-1)*SL(1);
	    else
		    ONE   = sqrt(FN*FN - FM*FM);
		    TWO   = sqrt((FN - 1.0)^2 - FM*FM)/ONE;
		    THREE = (2.0*FN - 1.0)/ONE;
		    I     = K - N;
		    J     = K - 2*N + 1;
		    P(K)  = (FN + 1.0)*(THREE*SLAT/FN*P(I) - TWO/(FN - 1.0)*P(J));
		    Q(K)  = THREE*(SLAT*Q(I) - CLAT/FN*P(I)) - TWO*Q(J);
        end
    end
    
    % Synthesize field components
    ONE   = (agh(L) + dgh(L)*T)*RR;
    if M == 0 %10,9,10
        X = X + ONE*Q(K);
	    Z = Z - ONE*P(K);
	    L = L + 1;
    else
	    TWO   = (agh(L+1) + dgh(L+1)*T)*RR;
	    THREE = ONE*CL(M) + TWO*SL(M);
	        X = X + THREE*Q(K);
	        Z = Z - THREE*P(K);
	    if CLAT > 0 %12,12,11
		    Y = Y+(ONE*SL(M)-TWO*CL(M))*FM*P(K)/((FN + 1.0)*CLAT);
	    else
		    Y = Y + (ONE*SL(M) - TWO*CL(M))*Q(K)*SLAT;
	    end
	        L = L + 2;
    end
	        M = M + 1;
end

%-------------------------------------------------------------------------
% Step 6: transform to geodetic if needed, compute total field
ONE   = X;
X     = X*CD +  Z*SD;
Z     = Z*CD - ONE*SD;
T     = sqrt(X*X + Y*Y + Z*Z);
J     = [X,Y,Z,T];
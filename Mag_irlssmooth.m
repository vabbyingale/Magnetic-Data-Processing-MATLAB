function y = Mag_irlssmooth(x,tau,ord)
    % irlssmooth: Iteratively Reweighted Least-Squares Smoother.
    %
    %   y = irlssmooth(x, tau) smooths vector x with a smoothing parameter tau,
    %   controlling the tradeoff between fidelity to data vs smoothness.
    %
    %   y = irlssmooth(x, tau, ord) uses an optional derivative order `ord`:
    %     - ord = 0..4 (default 2)
    %     - higher ord => sharper frequency knee, less damping
    %
    % For a given ORD, the unaffected space of TAU and L depends
    % on the machine and MATLAB version, so the rightmost constants in the
    % table below are approximate.
    %
    %                                      Approximate
    %           Order                  Domain of L and TAU
    %           -----               -------------------------
    %             0                      L^2 * TAU  <  1e28
    %             1                  L^(2/3) * TAU  <  2.7e9
    %             2 (default)        L^(2/5) * TAU  <  5.3e5
    %             3                  L^(2/7) * TAU  <  13000
    %             4                  L^(2/9) * TAU  <  1800
    %
    %    Note that as TAU surpasses the input length, the output approaches a
    %    polynomial fit of order ORD.
    %
    %    See also LSSMOOTH
    %
    %
    % Written by James S. Montanaro, February 2015

%-------------------------------------------------------------------------
% Step 1: Internal irls control parameters

maxcount = 25;               % max # iterations
tol      = 1e-6;             % relative output change tolerance
geb      = 1.5;              % Gaussian error breakpoint, in sigma +/-
g        = 1.3;              % loop gain

% Derived parameter
p2 = (geb/.6745)^2;          % squared breakpoint relative to median

%-------------------------------------------------------------------------
% Step 2: Check inputs

if nargin < 2
    error('Not enough inputs. Usage: y = irlssmooth(x, tau, ord)')
elseif tau < 0
    error('TAU must not be negative.')
elseif nargin < 3
    ord=2;
elseif ord < 0 || ord > 4
    error('ORD must be in the range 0 to 4.')
else
    ord=round(ord);
end

%-------------------------------------------------------------------------
% Step 3: Check dimensions 
[m,n] = size(x);
if min(m, n) ~= 1
    error('Input sequence must be a vector.')
end
row = (m == 1);
if row
    x = x.';
    m = n;
end
if m < ord+1
    error(['At least ' int2str(ord+1) ' input points required.'])
end
ordstart = ord;

%-------------------------------------------------------------------------
% Step 4: Algorithm

tau0 = [4.000 3.416 3.404 3.411 3.417];     % Tau normalizers for orders 0:4

% Prepare to reduce ORD if calculation is ill-conditioned
warnID = 'MATLAB:rankDeficientMatrix';      % MATLAB warning ID
s0     =  warning('off',warnID);            % disable warning message

% If necessary, repeat calculations until conditioning is good
ill = true;                                 % pretend ill-conditioned to start loop
ord = ord+1;                                % increment to cancel 1st decrement

while ill
    ord=ord-1;                              % reduce order
    lastwarn('');                           % clear warning indicator
    
    % Compute differencing coefficients
    h = [-1 1];                             % 1st-order diff
    for i = 1:ord
        h = conv(h,[-1 1]);                 % higher-order diff
    end

    k = ord+1;                              % for convenience
    h = repmat(h,m-k,1);                    % diagonals-to-be
    wd = (tau/tau0(k))^k;                   % weight of differencing part
   
    % Compute initial LS solution
    w=ones(m,1);                            % uniform error weighting
    Ad=wd*spdiags(h,0:k,m-k,m);             % differencing part of matrix
    A=[spdiags(w,0,m,m); Ad];               % sparse matrix
    v=[w.*x; zeros(m-k,1)];                 % target vector
    y=A\v;                                  % LS solution
    [lastmsg,lastID]=lastwarn;              % get last warning ID
    ill=strcmp(lastID,warnID);              % ill if lastID = warnID
    
    % IRLS loop
    count=0;                                % initialize loop counter
    chg=10*std(y-x)/std(detrend(y));        % initial output change (inflated)
    while ~ill && chg > tol && count < maxcount

        count=count+1;                      % increment counter
        oldy=y;                             % remember old solution
        oldchg=chg;                         % remember old output change
        e2=abs(x-y).^2;                     % squared "error" vector
        a2=p2*median(e2);                   % squared error breakpoint
        neww=sqrt(a2./(a2+e2));             % new error weighting, pre-gain
        w=w.^(1-g).*neww.^g;                % new error weighting
        w=w/mean(w);                        % preserve ratio of mean(w) to wd
        A=[spdiags(w,0,m,m); Ad];           % new sparse matrix
        v=[w.*x; zeros(m-k,1)];             % new target vector
        y=A\v;                              % new LS solution
        [lastmsg,lastID]=lastwarn;          % get last warning ID
        chg=std(y-oldy)/std(detrend(y));    % new output change
        ill=strcmp(lastID,warnID) || chg >= oldchg;  % illness criteria
        %disp(chg)
    end
    %count
    
end
warning(s0);                                % restore warning state

% -------------------------------- Notices --------------------------------

if ord < ordstart                     % if order was reduced
    %disp(['Notice: IRLSSMOOTH order ORD was reduced from ' ...
    %      int2str(ordstart) ' to ' int2str(ord) ...
    %      ' for conditioning.'])
    %disp('  ')
end
if count == maxcount                  % if max count was reached
    warning(['Maximum number of iterations reached.  ' ...
             'Per-iteration relative output change = ' num2str(chg)])
end

% Convert back to row vector if necessary
if row
    y=y.';
end

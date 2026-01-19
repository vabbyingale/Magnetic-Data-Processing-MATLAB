function [Fit_Matrix] = Mag_fit_spin_comp(heading,n_comp,e_comp,d_comp, sentry_number)
    % A function that fits robust sine models to magnetic field components for spin correction 
    %
    % Syntax:
    %   Fit_Matrix = fit_spin_comp(heading, n_comp, e_comp, d_comp);
    %
    % Inputs:
    %   heading = Heading angles (degrees, typically 0 to 360)
    %   n_comp  = North component of magnetic field (nT)
    %   e_comp  = East component of magnetic field (nT)
    %   d_comp  = Down component of magnetic field (nT)
    %
    % Outputs:
    %   Fit_Matrix = 3x3 matrix where each row corresponds to [a, b, c]
    %                parameters from fit: a + b*sin(c + heading)
    %                for North, East, and Down respectively.

%-------------------------------------------------------------------------
% Step 1: initial guess and preparing data

% Fit guess
guessa=100;                                 % offset guess
guessb=500;                                 % amplitude guess
guessc=20;                                  % phase shift guess
guessd=500;
guesse=20;

% Remove means
n_comp = n_comp - mean(n_comp);             % detrend N
e_comp = e_comp - mean(e_comp);             % detrend E
d_comp = d_comp - mean(d_comp);             % detrend D

%-------------------------------------------------------------------------
% Step 2: create different fitting equations

f = fittype('a+b*sind(c+x)+d*sind(e+2*x)');             

% fit 1: standard least square fit
[~,~,fitinfo_n] = fit(heading,n_comp,f,'StartPoint',[guessa guessb guessc guessd guesse]);
[~,~,fitinfo_e] = fit(heading,e_comp,f,'StartPoint',[guessa guessb guessc guessd guesse]);
[~,~,fitinfo_d] = fit(heading,d_comp,f,'StartPoint',[guessa guessb guessc guessd guesse]);


% fit 2: exclude outliers based on residuals (>2 std)
% Get outliers from fitinfo struct
residuals_n = fitinfo_n.residuals;          
residuals_e = fitinfo_e.residuals;
residuals_d = fitinfo_d.residuals;

% Find indices for residuals
I_n = abs(residuals_n) > 2 * std(residuals_n);
I_e = abs(residuals_e) > 2 * std(residuals_e);
I_d = abs(residuals_d) > 2 * std(residuals_d);

% Build exclusion mask
outliers_n = excludedata(heading,n_comp,'indices',I_n);
outliers_e = excludedata(heading,e_comp,'indices',I_e);
outliers_d = excludedata(heading,d_comp,'indices',I_d);

% Refit exluding outliers
fit2n = fit(heading,n_comp,f,'StartPoint',[guessa guessb guessc guessd guesse],'Exclude',outliers_n);
fit2e = fit(heading,e_comp,f,'StartPoint',[guessa guessb guessc guessd guesse],'Exclude',outliers_e);
fit2d = fit(heading,d_comp,f,'StartPoint',[guessa guessb guessc guessd guesse],'Exclude',outliers_d);

% fit 3: Robust fits (automatically downweights outliers)
fit3n = fit(heading,n_comp,f,'StartPoint',[guessa guessb guessc guessd guesse],'Robust','on');
fit3e = fit(heading,e_comp,f,'StartPoint',[guessa guessb guessc guessd guesse],'Robust','on');
fit3d = fit(heading,d_comp,f,'StartPoint',[guessa guessb guessc guessd guesse],'Robust','on');

%-------------------------------------------------------------------------
% Step 3: Plotting

fig1 = figure('Position', [200, 200, 1200, 600]);
s1 = subplot(1,3,1);
plot(heading,n_comp,'k.','MarkerSize',2); hold on
h3=plot(fit3n,'m-');
set(h3,'lineWidth',2);
title('North'); xlabel('Heading (°)'); ylabel('Magnetic Field (nT)')
xlim([0 360])
xticks([0 90 180 270 360])
xticklabels({'0','90','180','270','360'})
ylim_val = max(abs(ylim));
ylim([-ylim_val ylim_val])
grid on
legend('Data','Robust sin curve','Location','southeast')
hold off

s2 = subplot(1,3,2);
plot(heading,e_comp,'k.','MarkerSize',2); hold on
h4=plot(fit3e,'m-');
set(h4,'lineWidth',2);
title('East'); xlabel('Heading (°)'); ylabel('')
xticks([0 90 180 270 360])
xticklabels({'0','90','180','270','360'})
ylim_val = max(abs(ylim));
ylim([-ylim_val ylim_val])
xlim([0 360])
grid on
legend('hide')
hold off

s3 = subplot(1,3,3);
plot(heading,d_comp,'k.','MarkerSize',2); hold on
h5=plot(fit3d,'m-');
set(h5,'lineWidth',2);
title('Down'); xlabel('Heading (°)'); ylabel('')
xlim([0 360])
xticks([0 90 180 270 360])
xticklabels({'0','90','180','270','360'})
ylim_val = max(abs(ylim));
ylim([-ylim_val ylim_val])
legend('hide')
grid on
hold off

splots = [s1 s2 s3];
for n = 1:length(splots)
    Mag_position_plot(splots(n))
end
Mag_save_plot_prompt(fig1, sprintf('01_S%s_spins_fit', sentry_number));

%-------------------------------------------------------------------------
% Step 4: Store fit parameters in matrix
% Each row is [a b c] for N, E, D fits

Fit_Matrix = [
    fit3n.a fit3n.b fit3n.c; 
    fit3e.a fit3e.b fit3e.c; 
    fit3d.a fit3d.b fit3d.c];
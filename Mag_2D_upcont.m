clear all; close all; clc;

% Load files
data_mag     = load('File_forward_anom.txt');
data_depth   = load('File_forward_depth.txt');
data_forward = load('File_forward_model.txt');

% create meshgrid from unique x, y values
x = unique(data_mag(:,1));
y = unique(data_mag(:,2));
[X, Y] = meshgrid(x, y);

nx = length(x(:));
ny = length(y(:));

% interpolate to nearest point for given meshgrid
Z_mag = griddata(data_mag(:,1), data_mag(:,2), data_mag(:,3), X, Y, 'nearest');
Z_depth = griddata(data_depth(:,1), data_depth(:,2), data_depth(:,3), X, Y, 'nearest');
Z_forward = griddata(data_forward(:,1), data_forward(:,2), data_forward(:,3), X, Y, 'nearest');

% grid spacing from given file
dx = abs(x(2) - x(1));
dy = abs(y(2) - y(1));
zlev = ceil(max(Z_depth(:)));       % next integer to max depth (fortran input)
pad_factor = 0.1;                   % to avoid boundary problem of 2D grid

% Estimate ws based on max depth and grid spacing

lambda_cutoff = max(Z_depth(:))+abs(min(Z_depth(:)));
ws_x = lambda_cutoff / dx;
ws_y = lambda_cutoff / dy;
ws_guess = 3*mean([ws_x, ws_y]);        % emperical calculation to the best guess

% Estimate wl based on dominant anomaly wavelength from 2D FFT 
F = abs(fft2(Z_mag));
F_shift = fftshift(F);

% compute corresponding wavenumbers
kx = (-nx/2:nx/2-1)/(nx*dx);
ky = (-ny/2:ny/2-1)/(ny*dy);
[kxx, kyy] = meshgrid(kx, ky);

% Find dominant wavelength
[~, idx_max] = max(F_shift(:));
[ky_idx, kx_idx] = ind2sub(size(F_shift), idx_max);
lambda_dom = 1 / sqrt(kxx(ky_idx,kx_idx)^2 + kyy(ky_idx,kx_idx)^2);
wl_guess = lambda_dom / mean([dx, dy]);   

fprintf('Estimated wl ~ %d, ws ~ %d \n', wl_guess, ws_guess);

% compute upward continuation with guess wl, ws
[m_lev] = Mag_guspi_upward(Z_mag, Z_depth, dx, dy, wl_guess, ws_guess, zlev, pad_factor);

% Residual analysis between upward continued and forward model
diff = Z_forward - m_lev;
diff_flat = diff(:);
mu = mean(diff_flat);
sigma = std(diff_flat);

% Plot the results
zmin = round(min(Z_mag(:)) / 100) * 100;
zmax = round(max(Z_mag(:)) / 100) * 100;
contourLevels = zmin:200:zmax;

figure('Position', [100, 100, 1400, 1200]) 
tiledlayout(2, 3, 'Padding', 'tight', 'TileSpacing', 'tight');

% Plot 1: Model depth
nexttile
ax1 = gca;
contourf(X, Y, flipud(Z_depth), 'LineStyle', 'none'); 
colormap(ax1, 'parula')
colorbar(ax1);
clim(ax1, [min(Z_depth(:)), max(Z_depth(:))]); 
axis equal tight
Mag_position_plot(ax1)
set(gca,'YDir','Reverse')
title('Depth (km)');

% Plot 2: Original field
nexttile
ax2 = gca;
contourf(X, Y, flipud(Z_mag), 50, 'LineStyle', 'none'); 
hold on
[C2, h2] = contour(X, Y, flipud(Z_mag), contourLevels, 'k', 'LineWidth', 1);
clabel(C2, h2, 'FontSize', 10, 'Color', 'k');
colormap(ax2, 'jet');
clim(ax2, [min(Z_mag(:)), max(Z_mag(:))]); 
axis equal tight
colorbar;
Mag_position_plot(gca);
set(gca,'YDir','Reverse')
title('Original Magnetic Field');

% Plot 3: Forward field
nexttile
ax3 = gca;
contourf(X, Y, flipud(Z_forward), 50, 'LineStyle', 'none'); 
hold on
[C3, h3] = contour(X, Y, flipud(Z_forward), contourLevels, 'k', 'LineWidth', 1);
clabel(C3, h3, 'FontSize', 10, 'Color', 'k');
colormap(ax3, 'jet');
clim(ax3, [min(Z_mag(:)), max(Z_mag(:))]); 
colorbar;
axis equal tight
Mag_position_plot(gca);
set(gca,'YDir','Reverse')
title('Forward Model');

% Plot 4: Histogram of residuals
nexttile
ax_hist = gca;
max_abs = ceil(max(abs(diff_flat)));
edges = -max_abs:1:max_abs;
histogram(diff_flat, 'BinEdges', edges, 'Normalization', 'count');
xlabel('Residual = forward - upcont (nT)')
ylabel('Counts')
title(sprintf('\\mu = %.3f, \\sigma = %.3f', mu, sigma))
xlim([-max_abs, max_abs])
grid on
Mag_position_plot(ax_hist)
axis(ax_hist, 'square')
ylim padded

% Plot 5: Spatial residuals
residual_grid = reshape(diff_flat, size(Z_forward));
nColors = 256;
r1 = linspace(0, 1, nColors/2)'; 
g1 = linspace(0, 1, nColors/2)'; 
b1 = linspace(1, 1, nColors/2)'; 
r2 = linspace(1, 1, nColors/2)';
g2 = linspace(1, 0, nColors/2)';
b2 = linspace(1, 0, nColors/2)';
colorMap = [r1, g1, b1; r2, g2, b2];

cn = max(abs([min(flipud(Z_forward - m_lev)), max(flipud(Z_forward - m_lev))]));

nexttile
ax2 = gca;
contourf(X, Y, flipud(Z_forward - m_lev), 50, 'LineStyle', 'none'); 
colormap(ax2, colorMap)
clim(ax2, [-cn cn])
colorbar(ax2);
axis equal tight
Mag_position_plot(ax2)
set(gca,'YDir','Reverse')
title('Difference: Forward - Upcont');

% Plot 6: Upward continuation
nexttile
ax4 = gca;
contourf(X, Y, flipud(m_lev), 50, 'LineStyle', 'none'); 
hold on
[C1, h1] = contour(X, Y, flipud(m_lev), contourLevels, 'k', 'LineWidth', 1);
clabel(C1, h1, 'FontSize', 10, 'Color', 'k');
colormap(ax4, 'jet');
clim(ax4, [min(Z_mag(:)), max(Z_mag(:))]); 
colorbar; 
axis equal tight
xlim([min(x), max(x)])
set(gca,'YDir','Reverse')
Mag_position_plot(gca);
title(sprintf('Upcont at zlev = %.1f (wl=%.2f, ws=%.2f)', zlev, wl_guess, ws_guess));
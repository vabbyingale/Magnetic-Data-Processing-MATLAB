function [Filtered] = Mag_plot_filter_components(time, x_comp, y_comp, z_comp, win_size, title1, title2, title3, type)
    % A function to visualize three-component time series data with optional filtering.
    %
    % Syntax:
    %   Filtered = plot_filter_components(time, x_comp, y_comp, z_comp, ...
    %               win_size, title1, title2, title3, type)
    %
    % Inputs:
    %   time     - time vector
    %   x_comp   - x component data (e.g., North)
    %   y_comp   - y component data (e.g., East)
    %   z_comp   - z component data (e.g., Down)
    %   win_size - window size for filters (in samples)
    %   title1   - title for first subplot
    %   title2   - title for second subplot
    %   title3   - title for third subplot
    %   type     - filtering type:
    %              0: No filter (just plot raw)
    %              1: Moving median filter
    %              2: Moving median + IRLS smoothing
    %
    % Outputs:
    %   Filtered - matrix of filtered components [X Y Z]

%-------------------------------------------------------------------------
% Step 1: moving median filter
if type == 1                            
    x_comp_m = movmedian(x_comp, win_size);
    y_comp_m = movmedian(y_comp, win_size);
    z_comp_m = movmedian(z_comp, win_size);
    Filtered = [x_comp_m y_comp_m z_comp_m];
    
    figure('Position', [200, 200, 1200, 700])
    t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); 
    
    % subplot 1
    nexttile 
    plot(time, x_comp,   'k-', 'LineWidth', 1.5); hold on; 
    plot(time, x_comp_m, 'm-', 'LineWidth', 1.0); hold on;
    xlim([min(time) max(time)])
    xticklabels([])
    title(title1); grid on
    %ylim([25000 28700])
    legend({'Raw', 'Median'}, 'Orientation', 'horizontal')
    Mag_position_plot(gca)
    
    % subplot 2
    nexttile
    plot(time, y_comp,   'k-', 'LineWidth', 1.5); hold on; 
    plot(time, y_comp_m, 'm-', 'LineWidth', 1.0); hold on;
    xlim([min(time) max(time)])
    xticklabels([]); grid on
    ylabel('Magnetic Field (nT)')
    title(title2)
    %ylim([2000 8600])
    Mag_position_plot(gca)
    
    % subplot 3
    nexttile
    plot(time, z_comp,   'k-', 'LineWidth', 1.5); hold on; 
    plot(time, z_comp_m, 'm-', 'LineWidth', 1.0); hold on;
    xlim([min(time) max(time)])
    title(title3); grid on
    %ylim([-16200 -10000])
    Mag_position_plot(gca)

%-------------------------------------------------------------------------
% Step 2: moving median + L2 norm smoothing filter
elseif type == 2   
    x_comp_m = movmedian(x_comp, win_size);
    y_comp_m = movmedian(y_comp, win_size);
    z_comp_m = movmedian(z_comp, win_size);
    x_comp_i  = Mag_irlssmooth(x_comp, win_size, 2);
    y_comp_i  = Mag_irlssmooth(y_comp, win_size, 2);
    z_comp_i  = Mag_irlssmooth(z_comp, win_size, 2);
    Filtered  = [x_comp_i y_comp_i z_comp_i];
    
    figure('Position', [200, 200, 1200, 700])
    t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); 

    % subplot 1
    nexttile 
    plot(time, x_comp,   'k-', 'LineWidth', 1.5); hold on;
    plot(time, x_comp_m, 'm-', 'LineWidth', 1.0); hold on;
    plot(time, x_comp_i, 'b-', 'LineWidth', 1.0);
    xlim([min(time) max(time)])
    xticklabels([])
    title(title1); grid on
    legend({'Raw', 'Median', 'L2'}, 'Orientation', 'horizontal')
    Mag_position_plot(gca)

    % subplot 2
    nexttile 
    plot(time, y_comp,   'k-', 'LineWidth', 1.5); hold on;
    plot(time, y_comp_m, 'm-', 'LineWidth', 1.0); hold on;
    plot(time, y_comp_i, 'b-', 'LineWidth', 1.0);
    xlim([min(time) max(time)])
    xticklabels([]); grid on
    ylabel('Magnetic Field (nT)')
    title(title2)
    Mag_position_plot(gca)

    % subplot 3
    nexttile 
    plot(time, z_comp,   'k-', 'LineWidth', 1.5); hold on;
    plot(time, z_comp_m, 'm-', 'LineWidth', 1.0); hold on;
    plot(time, z_comp_i, 'b-', 'LineWidth', 1.0);
    xlim([min(time) max(time)])
    title(title3); grid on
    Mag_position_plot(gca)

%-------------------------------------------------------------------------
% Step 3: no filter
elseif type == 0                         
    figure('Position', [200, 200, 1200, 700])
    t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact'); 

    % subplot 1
    nexttile 
    plot(time,x_comp,'k-','LineWidth',0.5); hold on;
    xlim([min(time) max(time)])
    title(title1)
    xticklabels([]); grid on
    legend({'Raw'},'Orientation','horizontal')
    Mag_position_plot(gca)

    % subplot 2
    nexttile
    plot(time,y_comp,'k-','LineWidth',0.5); hold on;
    xlim([min(time) max(time)])
    ylabel('Magnetic Field (nT)')
    title(title2)
    xticklabels([]); grid on
    Mag_position_plot(gca)

    % subplot 3
    nexttile
    plot(time,z_comp,'k-','LineWidth',0.5); hold on;
    xlim([min(time) max(time)])
    title(title3); grid on
    Mag_position_plot(gca)
end
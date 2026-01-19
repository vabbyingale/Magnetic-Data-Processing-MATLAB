function Mag_plot_upcont_comparison(sentry_number)
    % Plot anomaly, upward continued fields, residuals, and depth for two methods.
    % User provides upward continuation level (in km).

    % -------------------------------------------------------------------------
    % Step 1: Get Sentry ID
    sentry_id = Mag_extract_sentry_id();

    % -------------------------------------------------------------------------
    % Step 2: Ask user for upward continuation level (e.g., 2.45 km)
    depth_km = input('Enter upward continuation depth (in km, e.g., 2.45): ');

    % Format depth for filenames (e.g., '2.45' → '2.45km')
    depth_str = sprintf('%.2fkm', depth_km);

    % -------------------------------------------------------------------------
    % Step 3: Construct filenames
    guspi_file = sprintf('S%s_guspi_upcont_%s.txt', sentry_number, depth_str);
    pk_file    = sprintf('S%s_PK_upcont_%s.txt', sentry_number, depth_str);

    % -------------------------------------------------------------------------
    % Step 4: Load data
    if ~isfile(guspi_file) || ~isfile(pk_file)
        error('One or both files do not exist:\n%s\n%s', guspi_file, pk_file);
    end

    guspi_data = load(guspi_file);
    pk_data    = load(pk_file);

    % Columns expected:
    % 1: distance (km)
    % 4: depth (m)
    % 5: observed anomaly
    % 6: upward continued anomaly

    % -------------------------------------------------------------------------
    % Step 5: Calculate residuals
    res_guspi = guspi_data(:,8) - guspi_data(:,12);
    res_pk    = pk_data(:,8)    - pk_data(:,12);

    % -------------------------------------------------------------------------
    % Step 6: Plot
    figure('Position',[200, 200, 1000, 900]);
    h = 0.268; gap = 0.038; startY = 0.705;

    % === First subplot: Anomalies ===
    s1 = subplot(3,1,1); set(s1, 'Position', [0.12, startY, 0.8, h])
    plot(guspi_data(:,1), guspi_data(:,8), 'LineWidth', 1.5, 'DisplayName', 'Original Anomaly'); hold on;
    plot(guspi_data(:,1), guspi_data(:,12), 'LineWidth', 1.5, 'DisplayName', 'Guspi UpCont'); hold on;
    plot(pk_data(:,1), pk_data(:,12), 'LineWidth', 1.5, 'DisplayName', 'PK UpCont');
    ylabel('Magnetic Anomaly (nT)');
    title('Comparing upward continued methods');
    legend; xticklabels([])
    grid on;

    % === Second subplot: Residuals ===
    s2 = subplot(3,1,2); set(s2, 'Position', [0.12, startY-(h+gap), 0.8, h])
    plot(guspi_data(:,1), res_guspi, 'Color', [128/255, 128/255, 128/255], 'LineWidth', 1.5, 'DisplayName', 'Guspi Residual'); hold on
    plot(pk_data(:,1), res_pk, 'Color', [50/255, 50/255, 50/255], 'LineWidth', 1.5, 'DisplayName', 'PK Residual'); hold on
    ylabel('Residual (nT)');
    legend; xticklabels([])
    grid on;

    % === Third subplot: Depth ===
    s3 = subplot(3,1,3); set(s3, 'Position', [0.12, startY-2*(h+gap), 0.8, h])
    plot(guspi_data(:,1), guspi_data(:,4), 'LineWidth',1.5, 'Color', [0/255, 0/255, 255/255], 'DisplayName','Depth'); hold on
    yline(depth_km, '--k', 'LineWidth', 1.2, 'DisplayName','Upcont Depth');
    set(gca, 'YDir','reverse');  % depth increases downward
    xlabel('Distance (km)'); legend
    ylabel('Depth (km)');
    ylim([depth_km-0.05 max(abs(guspi_data(:,4)))+0.05])
    grid on;

    splots = [s1 s2 s3];
    for n = 1:length(splots)
        Mag_position_plot(splots(n))
    end

end

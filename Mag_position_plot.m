function Mag_position_plot(ax)
    % A function for better visualization of plots
    
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
ax.XAxis.TickLength = [0.005 0.005];  
ax.YAxis.TickLength = [0.005 0.005]; 
ax.TickDir = 'out';  % Set ticks to point outward
ax.FontSize = 16; ax.FontName = 'Helvetica';
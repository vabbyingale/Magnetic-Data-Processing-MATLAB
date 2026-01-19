function Mag_save_plot_prompt(figHandle, defaultName)
%SAVE_PLOT_PROMPT Ask user if they want to save the current plot as PNG.
%
% figHandle   - Handle to the figure
% defaultName - Default filename (without extension)

    drawnow;  % Ensure plot is rendered
    choice = input('Save this plot as PDF? (y/n): ', 's');
    if strcmpi(choice, 'y')
        if nargin < 2 || isempty(defaultName)
            defaultName = sprintf('plot_%s', datestr(now, 'yyyymmdd_HHMMSS'));
        end
        filename = sprintf('%s.pdf', defaultName);
        saveas(figHandle, filename);
        fprintf('Plot saved as: %s\n', filename);
        fprintf('\n')
    end
    close(figHandle);
end

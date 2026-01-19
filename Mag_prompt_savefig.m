function Mag_prompt_savefig(fig_handle, fig_name, sentry_number, output_dir)
    % Ask user whether to save the figure. If yes, saves to specified directory.
    %
    % Inputs:
    %   fig_handle     - handle to the figure to save
    %   fig_name       - base filename (without extension)
    %   sentry_number  - e.g., '744' (will prepend 'S744_' to filename)
    %   output_dir     - optional: directory to save into (default: current folder)

if nargin < 4
    output_dir = '.';  % default to current folder
end

% Create prefix
prefix = sprintf('Fig_S%s_', sentry_number);

% Full figure filename
full_fig_name = [prefix fig_name];

% Ask user
choice = input(sprintf('\nSave figure "%s"? (y/n): ', full_fig_name), 's');

if strcmpi(choice, 'y')
    % Ensure output folder exists
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Save as both .png and .fig
    saveas(fig_handle, fullfile(output_dir, [full_fig_name, '.png']));

    fprintf('Saved: %s.png\n', full_fig_name);
else
    fprintf('Skipped saving figure "%s".\n', full_fig_name);
end


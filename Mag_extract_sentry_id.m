function sentry_id_str = Mag_extract_sentry_id()
    % Extract sentry dive ID from local file name
    %
    % Syntax:
    %   sentry_id_str = extract_sentry_id()
    % 
    % Output:
    %   sentry_id_str - string with the 3-digit Sentry number (e.g., '744'),
    %   or 'XXX' if not found.

% -------------------------------------------------------------------------
% Step 1: Search for files matching naming convention
d = dir('Seg_sentry*.txt');

% -------------------------------------------------------------------------
% Step 2: Attempt to extract Sentry ID using regex
if ~isempty(d)
    tokens = regexp(d(1).name, 'Seg_sentry(\d{3})', 'tokens', 'once');
    if ~isempty(tokens)
        sentry_id_str = tokens{1};
        return;
    end
end

% -------------------------------------------------------------------------
% Step 3: Fallback value if no match found
sentry_id_str = 'XXX';  % fallback


function flag = Mag_ask_yesno(prompt)
    % A function to prompt user for yes/no response
    %
    % Syntax:
    %   flag = ask_yesno(prompt) 
    %
    % Input:
    %   prompt - A string that will be displayed to the user as a prompt.
    %
    % Output:
    %   flag   - A logical value: true (1) if input is 'y' or 'Y', false (0) otherwise.

reply = input(prompt, 's');
flag = strcmpi(reply, 'y');
% Set a passcode either directly by supplying the passcode matrix otherwise
% a file choose will open for you.
function passcode=setPasscode(passcode)    
    % If left unspecified, ask user to supply a movie clip in '.mov' format
    if ~exist('passcode', 'var')
        [filename, pathname, ~] = uigetfile('*.mov', 'Load New Password Video');
        if isequal(filename,0) || isequal(pathname,0)
            passcode = [];
            return;
        end
        filename = [pathname filename];
        
        % == interpret video of gesture to passcode ==
        passcode=interpret_gesture(filename);
    end
    
    [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
    passcodefilename = [curr_path filesep 'currpasscode.mat'];
    save(passcodefilename, 'passcode'); 
end
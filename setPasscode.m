function passcode=setPasscode(passcode)
    passcode = [];
    
    % If left unspecified, ask user to supply a movie clip in '.mov' format
    if ~exist('passcode', 'var')
        [filename, pathname, ~] = uigetfile('*.mov', 'Load New Password Video');
        if isequal(filename,0) || isequal(pathname,0)
            return;
        end
        filename = [pathname filename];
        passcode=interpret_gesture(filename);
    end
    
    [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
    passcodefilename = [curr_path filesep 'currpasscode.mat'];
    save(passcodefilename, 'passcode'); 
end
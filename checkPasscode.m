% Check a passcode by supplying the passcode video file or otherwise
% a file chooser will open for you.
function passcode=checkPasscode(filename)
    passcode = [];
    
    % If left unspecified, ask user to supply a movie clip in '.mov' format
    if ~exist('filename', 'var')
        [filename, pathname, ~] = uigetfile('*.mov', 'Load New Password Video');
        if isequal(filename,0) || isequal(pathname,0)
            return;
        end
        filename = [pathname filename];
        
        
        % == interpret video of gesture to passcode ==
        passcode=interpret_gesture(filename);
    end
    
    [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
    passcodefilename = [curr_path filesep 'currpasscode.mat'];
    file=load(passcodefilename); 
    real_passcode = file.passcode;
    
    if any(size(passcode) ~= size(real_passcode)) || (any(any(real_passcode~=passcode)))
        msgbox('Wrong passcode.');
        passcode
        real_passcode
    else
        msgbox('Success!');
    end
end
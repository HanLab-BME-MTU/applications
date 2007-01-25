function success = startupTomlab
%STARTUPTOMLAB is a utility function that starts up tomlab if necessary
%
% SYNOPSIS: startupTomlab
%
% INPUT none
%
% OUTPUT success : true if everything's fine
%
% REMARKS
%
% created with MATLAB ver.: 7.3.0.298 (R2006b) on 
%
% created by: jdorn
% DATE: 25-Jan-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = false; % init success. Try not to print output argument

% Tomlab can only be run from linux
if ispc
    disp('To run tomlab from the server, please switch to linux')
else

    % check for 'tomRun' in the path
    trPath = which('tomRun');

    % if tomRun is somewhere, and if it is in a path that contains "tomlab",
    % we're fine. Otherwise, startup tomlab
    if ~isempty(trPath) && any(findstr(lower(trPath),'tomlab'))
        % all fine
        s = true;
    else
        try
            oldDir = cd('/opt/tomlab');
            startup;
            cd(oldDir);
            s = true;
        catch
            % if error, it's bad
            disp(sprintf('Error starting tomlab: %s',lasterr))
        end
    end
end

% report success
if nargout > 0
    success = s;
end
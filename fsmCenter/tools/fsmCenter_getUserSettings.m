function userDir=fsmCenter_getUserSettings
% fsmCenter_getUserSettings guides the user to selet a directory where to store its personal files
%
% SYNOPSIS      userDir=fsmCenter_getUserSettings
%
% INPUT         none
%
% OUTPUT        userDir : current user-sepcified directory; if the function fails to recover
%                         the directory, userDir=[] is returned
%
% DEPENDENCES   fsmMain uses { }
%               fsmMain is used by { fsmCenter }
%
% Aaron Ponti, March 11th, 2004

% Initialize userDir
userDir=[];

% Look for HOME directory
homeDir=getenv('HOME');
if isempty(homeDir)
    return
end

% Check whether a settings file already exists
if ispc==1
    iniFileName=[homeDir,filesep,'fsmWin.ini'];
else
    iniFileName=[homeDir,filesep,'fsmUnix.ini'];
end    
if exist(iniFileName)==2
    [label,currentDir]=textread(iniFileName,'%s %s');
    currentDir=char(currentDir); % Change of type
    if isempty(label) | strcmp(label,'USERDIR')==0
        % Failure
        return
    else
        if exist(currentDir)==7
            % Success
            userDir=currentDir;
        else
            % Failure
            return
        end
    end
else
    % Failure - return
    return
end

function status=fsmToolsPrepareWorkDir(userPath)
% fsmToolsPrepareWorkDir creates the work directory specified by fsmParam
%
% SYNOPSIS      status=fsmToolsPrepareWorkDir(userPath)
%
% INPUT         userPath : path to be generated
%                          In 'path', following subdirectories are created
%                                       bwMask
%                                       cands
%                                       gapList
%                                       kinScore
%                                       locMax
%                                       movies
%                                       vectors
%
% OUTPUT        status   : returns 1 if the path has been generated correctly
%
% DEPENDENCES   
%               
%
% Aaron Ponti, October 22nd, 2002

% Set initial status to 0
status=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PREPARE WORKING DIRECTORY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store current directory
oldDir=cd;

% Default
st=1;

% Create ROOT directory
if exist(userPath)~=7
    % Directory does not exist - create it
    if userPath(2)==':';   % Windows
        % Drive letter specified
        st=mkdir(userPath(1:3),userPath(4:end));
    else
        st=mkdir(userPath);
    end
end

% Check status
if st~=1
    return;
end
    
% Change to current directory
cd(userPath);

% Add subdirectories
if exist([userPath,filesep,'bwMask'])~=7
    st=mkdir('bwMask');
    % Check status
    if st~=1
        return;
    end
end
if exist([userPath,filesep,'cands'])~=7
    st=mkdir('cands');
    % Check status
    if st~=1
        return;
    end
end
if exist([userPath,filesep,'gapList'])~=7
    st=mkdir('gapList');
    % Check status
    if st~=1
        return;
    end
end
if exist([userPath,filesep,'kinScore'])~=7
    st=mkdir('kinScore');
    % Check status
    if st~=1
        return;
    end
end
if exist([userPath,filesep,'locMax'])~=7
    st=mkdir('locMax');
    % Check status
    if st~=1
        return;
    end
end
if exist([userPath,filesep,'movies'])~=7
    st=mkdir('movies');
    % Check status
    if st~=1
        return;
    end
end
if exist([userPath,filesep,'vectors'])~=7
    st=mkdir('vectors');
    % Check status
    if st~=1
        return;
    end
end
if exist([userPath,filesep,'tftraw'])~=7
    st=mkdir('tftraw');
    % Check status
    if st~=1
        return;
    end
end


% Set status to 1
status=1;
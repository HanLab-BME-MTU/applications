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

% Default
st=1;

% Create ROOT directory
if ~exist(userPath, 'dir')
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
if ~exist([userPath,filesep,'bwMask'], 'dir')
    st=mkdir('bwMask');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'cands'], 'dir')
    st=mkdir('cands');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'gapList'], 'dir')
    st=mkdir('gapList');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'kinScore'], 'dir')
    st=mkdir('kinScore');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'locMax'], 'dir')
    st=mkdir('locMax');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'movies'], 'dir')
    st=mkdir('movies');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'vectors'], 'dir')
    st=mkdir('vectors');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'links'], 'dir')
    st=mkdir('links');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'flow'], 'dir')
    st=mkdir('flow');
    % Check status
    if st~=1
        return;
    end
end
if ~exist([userPath,filesep,'pointFiles'], 'dir')
    st=mkdir('pointFiles');
    % Check status
    if st~=1
        return;
    end
end
% Set status to 1
status=1;
function mainDir=cdSimdata(option)
%change to simdata/home if not in a simdata-dir/home-dir already
%
%SYNOPSIS   mainDir=cdSimdata(option)
%
%INPUT      option: 0 always cd simdata
%                   1 if in a simdata directory: don't do anything
%                   2 if in a simdata directory and if there are sim*-files: move up one level
%
%OUTPUT     SIMDATA - directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mainDir=getenv('SIMDATA');
if isempty(mainDir)
    disp('no SIMDATA environment variable - going HOME instead')
    mainDir=[getenv('HOME'),filesep,'matlab'];
    if isempty(mainDir)
        h=errordlg('Set environment variables HOME and SIMDATA for your matlab home directory and the directory where you store the data, respectively');
        uiwait(h)
        return
    end
end

oldDir=pwd;

%if there is no option: switch to simdata-main if not already in some
%simdata-subdirectory
if nargin==0|isempty(option)
    option=1;
end

switch option
    case 0
        cd(mainDir);
    case 1
        if isempty(findstr(lower(mainDir),lower(oldDir))) %else we're in a good directory...
            cd(mainDir);
        elseif length(mainDir)>length(oldDir) %...or we are in a parent
            cd(mainDir);
        end 
    case 2
        if isempty(findstr(lower(mainDir),lower(oldDir)))
            cd(mainDir);
        else
            if ~isempty(dir('sim*'))
                cd .. %if in a project dir: move up one file to allow easier switching between projects
            end
        end 
end
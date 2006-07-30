function [mainDir,oldDir]=cdBiodata(option)
%change to biodata/home if not in a biodata-dir/home-dir already
%
%SYNOPSIS  [mainDir,oldDir]=cdBiodata(option)
%
%INPUT      option: 0 always cd biodata
%                   1 if in a biodata directory: don't do anything
%                   2 if in a biodata directory and if there are r3d-files: move up one level
%                   3 as 2, but if not in biodata directory: do not change to it
%                   4 do not change dir at all, just return mainDir
%                   5 goto tmpData/misteli (only works in LCCB)
%                   6 goto tmpData/jason (only works in LCCB)
%
%OUTPUT     mainDir BIODATA - directory
%           oldDir  previous directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mainDir=getenv('BIODATA');
if isempty(mainDir) || ~isdir(mainDir)
    mainDir=getenv('HOME');
    if isempty(mainDir) || ~isdir(mainDir)
        h=errordlg('Set environment variables HOME and BIODATA for your matlab home directory and the directory where you store the movies, respectively');
        uiwait(h)
        return
    end
end

oldDir=pwd;

%if there is no option: switch to biodata-main if not already in some
%biodata-subdirectory
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
            if ~isempty(dir('*.r3*'))
                cd .. %if in a project dir: move up one file to allow easier switching between projects
            end
        end 
    case 3
        if ~isempty(dir('*.r3*'))
            cd .. %if in a project dir: move up one file to allow easier switching between projects
        end
    case 4
        % do not do anything - just return mainDir
    case 5
        % go to tmpData/misteli
        
       
        
        % cd Biodata and then misteli
        cdBiodata(0);
        try
        cd ../tmpData/misteli
        catch
        end
        mainDir = pwd;
        
       
    case 6
        % go to tmpData/jason
        
        % cd Biodata
        cdBiodata(0);
        try
        cd ../tmpData/jason       
        catch
            % we're on the laptop
            cd c:\tmp
        end
        mainDir = pwd;
        
    case 7
        % go to tmpData/misteli, but stay if there already
        
        % if we're already there, go maybe one dir up
        if any(findstr(oldDir,'misteli'))
            % check whether we're in a movieDir
            if isempty(dir('*.STK'))
                % we're fine
                mainDir = pwd;
            else
                cd ..
                mainDir = pwd;
            end
        else
            cdBiodata(0);
        cd ../tmpData/misteli
        mainDir = pwd;
        end
end
function fsmCenter_setUserSettings
% fsmCenter_setUserSettings guides the user to selet a directory where to store its personal files
%
% SYNOPSIS      fsmCenter_setUserSettings
%
% INPUT         none
%
% OUTPUT        none
%
% DEPENDENCES   fsmMain uses { }
%               fsmMain is used by { fsmCenter }
%
% Aaron Ponti, March 11th, 2004

% Look for HOME directory
homeDir=getenv('HOME');
if isempty(homeDir)
    uiwait(msgbox('Please make sure you have set an ENV variable called HOME. User directory not set.','Error','modal'));
    return
end

% Flag
selectDir=0;

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
        uiwait(msgbox('Invalid .ini file. Please select another directory.','Error','modal'));
        selectDir=1;
    else
        if exist(currentDir)==7
            string=['Your current directory is [ ',currentDir,' ]. Do you want to change it?'];
            choice=questdlg(string,'User input requested','Yes','No','No');
            switch choice,
                case 'Yes', selectDir=1;
                case 'No', selectDir=0;
            end % switch
        else
            string=['Current directory [ ',currentDir,' ] does not exist. Please select another one.'];
            uiwait(msgbox(string,'Error','modal'));
            selectDir=1;
        end
    end
else
    % ini file does not exist. The user MUST select a directory.
    selectDir=1;
end    

if selectDir==1
    userDir=uigetdir('','Pick a directory where to store your settings...');
    if userDir==0
        uiwait(msgbox('User directory not set.','Warning','modal'));
        return
    end
    
    % Save
    fid=fopen(iniFileName,'w+');
    if fid==-1
        uiwait(msgbox('Cannot create file .ini in HOME. User directory not set.','Error','modal'));
        return
    end
    fprintf(fid,'%s %s','USERDIR',userDir);
    fclose(fid); 
    
    % Copy file
    replaceFile=1;
  
    % Path of the default fsmExpParams.txt file
    pathOfFsmMain=which('fsmMain.m');
    % Get path for fsmParam
    indx=find(pathOfFsmMain==filesep);
    indx=indx(length(indx));
    source=[pathOfFsmMain(1:indx),'fsmExpParams.txt'];
    
    % Path of the potential destination file
    destFileName=[userDir,filesep,'fsmExpParams.txt'];
    
    if strcmp(source,destFileName)==1
        % The user selected the main fsm directory, where the default fsmExpParams.txt
        % file is stored.
        uiwait(msgbox('You selected the default [ fsmExpParams.txt ] file. This is acceptable, but it is recommended that you select another location.','Info','modal'));
        return
    end

    if exist(destFileName)==2
        string=['A file [ fsmExpParams.txt ] already exists in the specified directory. Do you want to overwrite it with the default file?'];
        choice=questdlg(string,'User input requested','Yes','No','No');
        switch choice,
            case 'Yes', replaceFile=1;
            case 'No', replaceFile=0;
        end % switch
    end
    if replaceFile==1
              
        if isempty(source)
            uiwait(msgbox('Could not find default experiment database! Get fsmExpParams.txt from the repository and restart SpeckTackle.','Error','modal'));
            return
        end
        success=copyfile(source,destFileName);
        if success==1
            uiwait(msgbox('File [ fsmExpParams.txt ] successfully copied.','Info','modal'));
        else
            uiwait(msgbox('Copy failed.','Error','modal'));
        end
    end
   
end 

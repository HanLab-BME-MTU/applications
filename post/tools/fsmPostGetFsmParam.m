function [fsmParam,status]=fsmPostGetFsmParam(projDir)

% Set status flag to zero (error)
status=0;

% Initialize output
fsmParam=[];

if isempty(projDir)
    % Select fsmParam.mat
    [fsmParamName,projDir] = uigetfile(...
        {'fsmParam.mat','fsmParam.mat'},...
        'Select fsmParam.mat');
    if projDir==0
        % The user pressed on Cancel
        return;
    end
else
    fsmParamName='fsmParam.mat';
    if exist([projDir,filesep,fsmParamName],'file')==0
        errordlg('No fsmParam.mat file could be found in the passed project directory.','Error','modal');
        return
    end    
end

if(isa(projDir,'char') & isa(fsmParamName,'char'))
    loadedFile=load([projDir,filesep,fsmParamName]);
    if strcmp(fieldnames(loadedFile),'fsmParam')
        fsmParam=loadedFile.fsmParam;
    else
        errordlg('The loaded file is not a valid ''fsmParam.mat'' file.','Error','modal');
        return
    end
else
    return 
end

% Set status to 1 to indicate success
status=1;
function fsmBatchJobs
% fsmBatchJobs runs several FSM calculation in a batch job
%
% Run this function from the fsmCenter user interface.
%
% SYNOPSIS fsmBatchJobs
%
% INPUT    NONE
%
% OUTPUT   NONE
%
% REMARK   fsmbatchJobs needs fsmParam.mat files to be saved into their corresponding work
%          directory. You can used fsmGuiMain to set up experiments and to save their
%          parameters as fsmParam.mat files.

% Store current directory
oldDir=cd;

% Structure containing all directories
userPath=struct('directory','','firstImageFile','');

% Load all fsmParam.mat files for the batch file
dName='default';
counter=0;
while ~isempty(dName)
    % Select fsmParam.mat
    [fName,dName] = uigetfile(...
        {'fsmParam.mat','fsmParam.mat'},...
        'Select fsmParam.mat - Click Cancel to stop adding jobs');
    if(isa(fName,'char') & isa(dName,'char'))
        counter=counter+1;
        userPath(counter).directory=dName;
    else
        break
    end

    % Select first image
    [fileName,dirName] = uigetfile(...
        {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
        '*.tif','TIF files (*.tif)'
        '*.tiff','TIFF files (*.tiff)'
        '*.jpg;','JPG files (*.jpg)'
        '*.jpeg;','JPEG files (*.jpeg)'
        '*.*','All Files (*.*)'},...
        'Select first image');
    if(isa(fName,'char') & isa(dirName,'char'))
        userPath(counter).firstImageFile=[dirName,fileName];
    else
        error('Image file not specified');
    end
end
 
% Start batch job
if isempty(userPath(1).directory)
    return
end

% Start gui
h=batchGUI;
set(h,'Name','Running...');
strH=get(h,'Children');
nExp=length(userPath);
set(strH(1),'String',num2str(nExp));
set(strH(4),'String','Running batch job');
set(strH(3),'String','0');
set(strH(2),'String','of');


for i=1:nExp
    
    % Update gui
    set(strH(3),'String',num2str(i));
    
    cd(userPath(i).directory);
    if ~exist('fsmParam.mat')
        error('No fsmParam.mat found in current directory');
    end
    load fsmParam.mat
    
    % Check fsmParam
    if isempty(fsmParam.specific.fileList)
        % This experiment has never been run before - fill the specific subfield
        fsmParam.specific.imgSize=size(imread(userPath(i).firstImageFile));
        % Number of image files
        fsmParam.specific.imageNumber=fsmParam.main.imgN;% Get filename list
        % Get image file list
        outFileList=getFileStackNames(userPath(i).firstImageFile);
        [path,body,no,ext]=getFilenameBody(char(outFileList(end)));
        % Prepare string number format
        s=length(num2str(no));
        strg=sprintf('%%.%dd',s);
        fsmParam.specific.formString=strg;
        % Store file list
        outFileList=outFileList(1:fsmParam.main.imgN);
        fsmParam.specific.fileList=char(outFileList);
        % Recover the index of the first and last image treated
        [path,body,no,ext]=getFilenameBody(char(outFileList(1)));
        fsmParam.specific.firstIndex=str2num(no);
        [path,body,no,ext]=getFilenameBody(char(outFileList(end)));
        fsmParam.specific.lastIndex=str2num(no);
        % Bleaching correction
        fsmParam.specific.intCorrFactors=ones(1,length(outFileList));
    end

    % Add field for batch jobs
    fsmParam.batchJob=1;
    
    % Perform calculation
    fsmMain(fsmParam); 
    
end

% Update gui
if nExp>1
    str=[num2str(nExp),' batch jobs run.'];
else
    str='1 batch job run.';
end 
set(strH(4),'String',str);
set(strH(3),'String','');
set(strH(2),'String','');
set(strH(1),'String','');
set(h,'Name','All done!');

% Change back to old dir
cd(oldDir);

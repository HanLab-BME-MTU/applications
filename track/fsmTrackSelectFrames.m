function [fsmParam,status]=fsmTrackSelectFrames(fsmParam)
% fsmTrackSelectFrames updates fsmParam in the case the user wants to track the sum of successive partial preprocessings
%
% SYNOPSIS   fsmParam=fsmTrackSelectFrames(fsmParam)
%
% INPUT      fsmParam:   general parameter structure
%
% OUTPUT     fsmParam:   modified (if needed) parameter structure
%            status  :   reports a failure
%
% DEPENDENCES   fsmTrackSelectFrames uses {}
%               fsmTrackSelectFrames is used by { fsmTrackMain }
%
% Aaron Ponti, January 7th, 2004

% Set initial status to failure
status=0;

% Check input parameter
if nargin~=1
    error('Input parameter fsmParam expected');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% READ NEEDED PARAMETERS FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

userPath=fsmParam.main.path;
origFileList=fsmParam.specific.fileList;
origFirstIndex=fsmParam.specific.firstIndex;
origLastIndex=fsmParam.specific.lastIndex;
xmin=fsmParam.main.normMin;      % Lower intensity bound for intensity normalization
xmax=fsmParam.main.normMax;      % Upper intensity bound for intensity normalization
tracker=fsmParam.track.tracker;
imageNumber=fsmParam.specific.imageNumber;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% READ CANDS LIST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd([userPath,filesep,'cands']);

% Look for the first cands###.mat file.
allFiles=dir;
i=1;
while ~strncmp('cands',allFiles(i).name,5)
    i=i+1;
end

% Read all cands###.mat file names
outFileList=getFileStackNames(allFiles(i).name);

% Read first and last indices
[path,body,firstIndex]=getFilenameBody(outFileList{1});
[path,body,lastIndex]=getFilenameBody(outFileList{end});

% Ask for user input - Setup dialog 
first=str2double(firstIndex);
last=str2double(lastIndex);
if tracker==3
    frames=2;
else
    frames=1;
end

% Check that there exist enough preprocessed frames to track at all
if length(outFileList)<(frames+1)
    if frames==1
        uiwait(msgbox('At least 2 preprocessed images are needed for tracking. Please run the PREPROCESSING module again.','Error','error','modal'));
        return % This returns an error to fsmMain (status=0)
    else
        uiwait(msgbox('At least 3 preprocessed images are needed for tracking. Please run the PREPROCESSING module again.','Error','error','modal'));
        return % This returns an error to fsmMain (status=0)
    end
end

% Select frame pairs to be processed
titleGUI=['Processed frames (first contiguous series): ',firstIndex,' to ',lastIndex];
[uFirst,uLast]=fsmTrackSelectFramesGUI(first,last,frames,titleGUI);

% Check whether the dialog was closed (_CloseRequestFcn)
if uFirst==-1
    return % This will return an error (status=0)
end

% Update fsmParam - uFirst and uLast are global variables returned by the dialog
if uFirst~=origFirstIndex || uLast~=origLastIndex

    % Create image file list
    [path,body,firstImageIndex,ext]=getFilenameBody(origFileList{1});
    % Replace numerical extension saved in the original image file list with that of
    % the first preprocessed image
    imageFileName=[path,filesep,body,firstIndex,ext]; 
    outImageFileList=getFileStackNames(imageFileName);
    outImageFileList=outImageFileList(uFirst-str2double(firstIndex)+1:uLast-str2double(firstIndex)+1);

    % Update fsmParam
    num=uLast-uFirst+1;
    fsmParam.main.imgN=num;
    fsmParam.specific.imageNumber=num;
    fsmParam.specific.fileList=outImageFileList;
    fsmParam.specific.firstIndex=uFirst;
    fsmParam.specific.lastIndex=uLast;
    
    % factors=fsmPrepIntCorrFactors(outFileList,num,[xmin xmax]);
    factors=ones(1,num);   % The intensity correction has been removed; yet, the structure has been maintained
    fsmParam.specific.intCorrFactors=factors;
end

% Set status to success
status=1;



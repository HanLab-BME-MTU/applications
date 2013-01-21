function segmentFAsTIRFStack(firstImageFile,saveDir,thresholdMethod,...
    methodValue,filterNoise,filterBackground,minSize,plotRes,firstMaskFile)
%segmentFAsTIRFStack segments FAs in a stack of images
%
%SYNOPSIS segmentFAsTIRFStack(firstImageFile,saveDir,thresholdMethod,...
%    methodValue,filterNoise,filterBackground,minSize,plotRes,firstMaskFile)
%
%INPUT  firstImageFile: Location and name of first image file.
%       saveDir       : Directory where masks are to be saved.
%       thresholdMethod: 
%                   'otsu' for Otsu.
%                   'rosin' for Rosin.
%                   'minmax' for first minimum after first maximum.
%                   'prct' to use a certain percentile.
%                   'user' to use a threshold input by user.
%                   Optional. Default: 'otsu'.
%       methodValue: Needed only if thresholdMethod = 'prct' or 'user'.
%                    If 'prct', then this is the percentile to use.
%                    Optional. Default: 90.
%                    If 'user', then this is the threshold value to use.
%                    Optional. Default: 0.9.
%       filterNoise: Either 0 to not filter noise or filter sigma > 0 to
%                    filter noise.
%                    Optional. Default: 1.
%       filterBackground: Either 0 to not filter background or filter sigma
%                         > 0 to filter background.
%                         Optional. Default: 10.
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%       firstMaskFile: Location and name of first mask to use for segmentation.
%                   Optional. If not provided, the whole image domain is segmented.
%
%OUTPUT 
%
%Khuloud Jaqaman January 2013

%% Output

%% Input

%first image file
if nargin < 1 || isempty(firstImageFile)
    [fName,dirName] = uigetfile('*.tif','specify first image in the stack - specify very first image');
else
    if iscell(firstImageFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstImageFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
    numFrames = length(outFileList);
    
    %read images
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);
    imageStack = NaN(isx,isy,numFrames);
    imageStack(:,:,1) = currentImage;
    for iFrame = 2 : numFrames
        imageStack(:,:,iFrame) = imread(outFileList{iFrame});
    end
    
else %else, exit
    disp('--segmentFAsTIRFStack: Bad file selection');
    return
end

%directory for saving
if nargin < 2 || isempty(saveDir)
    saveDir = uigetdir([],'Choose directory to save results');
end

%thresholding method
if nargin < 3 || isempty(thresholdMethod)
    thresholdMethod = 'Otsu';
end

%method value
if nargin < 4 || isempty(methodValue)
    switch thresholdMethod
        case 'prct'
            methodValue = 90;
        case 'user'
            methodValue = 0.9;
        otherwise
            methodValue = [];
    end
end

%noise filtering
if nargin < 5 || isempty(filterNoise)
    filterNoise = 1;
end

%background filtering
if nargin < 6 || isempty(filterBackground)
    filterBackground = 10;
end

%minimum blob size
if nargin < 7 || isempty(minSize)
    minSize = 20;
end

%plot results
if nargin < 8 || isempty(plotRes)
    plotRes = 0;
end

%first mask file
if nargin < 9 || isempty(firstMaskFile)
    [fName,dirName] = uigetfile('*.tif','specify first mask in the stack - if no mask, press ''Cancel''');
else
    if iscell(firstMaskFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstMaskFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstMaskFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstMaskFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%read masks
maskStack = true(isx,isy,numFrames); %default is mask=1, i.e. segment whole image
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileListM = getFileStackNames([dirName,fName]);
    
    %read images
    for iFrame = 1 : numFrames
        maskStack(:,:,iFrame) = imread(outFileListM{iFrame});
    end
    
end

if ~logical(maskStack)
    error('Mask must be a logical image.');
end
    
%% Segmentation

%go to directory where results will be saved
cd(saveDir)

%save segmentation parameters
save('segmentParam','thresholdMethod','methodValue','filterNoise','filterBackground','minSize');

%make directory to save masks
mkdir masks
saveDir = [saveDir filesep 'masks'];

%go over all frames ...
for iFrame = 1 : numFrames
    
    %call segmentatin code
    maskBlobs = segmentFAsTIRF(imageStack(:,:,iFrame),thresholdMethod,...
        methodValue,filterNoise,filterBackground,minSize,plotRes,maskStack(:,:,iFrame));
    
    %save masks
    [~,fname,fno,fext]=getFilenameBody(outFileList{iFrame});
    filename = [saveDir filesep 'maskFA_' fname fno fext];
    imwrite(maskBlobs,filename,'tif');
    
end

%% ~~~ the end ~~~


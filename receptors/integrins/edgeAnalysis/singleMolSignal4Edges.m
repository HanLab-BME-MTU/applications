function singleMolSignal4Edges(firstEdgeFile,firstSMFile,dir2save,...
    stepSize,numSMFrames)
%singleMolSignal4Edges derives images from single molecule signal to be used in edge detection
%
%SYNOPSIS singleMolSignal4Edges(firstEdgeFile,firstSMFile,dir2save,...
%    stepSize,numSMFrames)
%
%INPUT  firstEdgeFile: Full path + file name of first edge image in stack.
%       firstSMFile  : Full path + file name of first single molecule
%                      image in stack.
%       dir2save     : Directory name where to save modified edge frames.
%       stepSize     : Number of single molecule frames between consecutive
%                      edge frames.
%       numSMFrames  : Number of single molecule frames after edge frame to
%                      use in generating new image
%
%OUTPUT A stack of compound images derived from single molecule data,
%       stored in the directory "dir2save".
%
%REMARKS Run this function after separating the edge frames from the single
%        molecule stack and reenumerating them.
%
%Khuloud Jaqaman, October 2012

%% Input

if nargin ~= 5
    error('singleMolSignal4Edges: Wrong number of input arguments');
end

%% Processing

%get information about edge files
[fpathE,fnameE,fnoE,fextE] = getFilenameBody(firstEdgeFile);
dirName = [fpathE filesep];
fName = [fnameE fnoE fextE];

%get all file names in stack
edgeFileList = getFileStackNames([dirName fName]);
numEdgeFiles = length(edgeFileList);

%get information about single molecule files list
[fpathSM,fnameSM,fnoSM,fextSM] = getFilenameBody(firstSMFile);

%get number of digits used for enumeration
numDigits = length(fnoSM);

formatString = ['%0' num2str(numDigits) 'i'];

%go over all edge images and derive "single molecule images"
for iEdge = 1 : numEdgeFiles
    
    %read edge image
    edgeImage = double(imread(edgeFileList{iEdge}));
    
    %decide which single molecule images to read
    smImageIndx = (iEdge-1)*stepSize + (2:numSMFrames);
    
    %go backwards if there are no single molecule images after edge frame
    numString = num2str(smImageIndx(1),formatString);
    imageName = [fpathSM filesep fnameSM numString fextSM];
    if ~exist(imageName,'file')
        smImageIndx = (iEdge-1)*stepSize + (-numSMFrames+2:0);
    end
    
    %read single molecule images
    smImage = repmat(zeros(size(edgeImage)),[1 1 numSMFrames-1]);
    for iSM = 1 : numSMFrames-1
        numString = num2str(smImageIndx(iSM),formatString);
        imageName = [fpathSM filesep fnameSM numString fextSM];
        smImage(:,:,iSM) = double(imread(imageName));
    end
    
    %put the single molecule images together
    smImage2Use = max(smImage,[],3);
    
    %save new single molecule image
    numString = num2str(iEdge,formatString);
    fileName = [dir2save filesep 'sm4edge_' fnameE numString fextE];
    imwrite(uint16(smImage2Use),fileName,'tif');
    
end

%% ~~~ the end ~~~
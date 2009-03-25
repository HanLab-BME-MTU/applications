function subRoiConfig
% SUBROICONFIG allows user to choose sub-regions of interest and get detection & tracking from original ROI


% ask user for project directory
roiDir=uigetdir(pwd,'Please select roi_x directory to use for making sub-rois');

% find images directory
homeDir=pwd;
cd(roiDir)
cd ..
imDir=[pwd filesep 'images'];

% load roiYX and roiMask
roiYX=load([roiDir filesep 'roiYX.mat']);
roiYX=roiYX.roiYX;
roiMask=imread([roiDir filesep 'roiMask.tif']);

% load detected features from feat directory
movieInfo=load([roiDir filesep 'feat' filesep 'movieInfo']);
movieInfoOld=movieInfo.movieInfo;

% load tracking data and convert tracksFinal to matrix
trackDir = [roiDir filesep 'track'];
[listOfFiles]=searchFiles('.mat',[],trackDir,0);
load([listOfFiles{1,2} filesep listOfFiles{1,1}])
[trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal);
clear trackedFeatureIndx

%get number of tracks and number of time points
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints=numTimePoints./8;
% get start and end columns for each frame
eC=8*[1:numTimePoints]';
sC=eC-7;

% extract coordinates of tracks
tracks.xCoord = trackedFeatureInfo(:,1:8:end); 
tracks.yCoord = trackedFeatureInfo(:,2:8:end);

% get list of images, read first image, and make RGB (gray)
[listOfImages]=searchFiles('.tif',[],imDir,0);
img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
[imL,imW]=size(img(:,:,1));
img=(img-min(img(:)))./(max(img(:))-min(img(:)));
img=repmat(img,[1,1,3]);


% make sub-roi directory under roi_x directory and go there
subRoiDir=[roiDir filesep 'subROIs'];
if isdir(subRoiDir)
    rmdir(subRoiDir,'s');
end
mkdir(subRoiDir);
cd(subRoiDir)

% input to varycolor is black (we don't use it since the image is dark;
% make it one more than the number of ROIs (max is 9)
cMap=varycolor(10);

% set cell boundary to white in composite image
[img2show]=addMaskInColor(img,roiMask,[1 1 1]);

roiCount=1; % counter for rois for current project
makeNewROI=1; % flag for making new roi

% store 1 in every pixel of sub_1, 2 in every pixel of sub_2, etc.
labelMatrix=zeros(size(roiMask));

% iterate til the user is finished
selectROI=1;
while makeNewROI==1 && roiCount<10
    % make new roi_n image/analysis directories
    currentRoiAnDir=[pwd filesep 'sub_' num2str(roiCount)];
    mkdir(currentRoiAnDir);

    tempRoi=[];
    while isempty(tempRoi)
        try
            % draw polygon to make mask
            [tempRoi,polyXcoord,polyYcoord]=roipoly(img2show);
        catch
            disp('Please try again.')
        end
    end
    close all

    
    % get intersection with max region in the cell
    tempRoi=tempRoi & roiMask;
    % shrink max region in the cell for next round by excluding current roi
    roiMask=roiMask-tempRoi;
    % fill in label matrix
    labelMatrix(tempRoi)=roiCount;
   
    % add the current roi to the composite image
    [img2show]=addMaskInColor(img2show,tempRoi,cMap(roiCount,:));
    
    % get coordinates of vertices (ie all pixels of polygon boundary)
    [y1,x1]=ind2sub([imL,imW],find(tempRoi,1)); % first pixel on boundary
    roiYX = bwtraceboundary(tempRoi,[y1,x1],'N'); % get all pixels on boundary
    % test to make sure roi can be reproduced ok - assume that new polygon
    % doesn't differ in area more than 10% of tempRoi
    [resultBW]=roipoly(imL,imW,roiYX(:,2),roiYX(:,1));
    if abs(sum(resultBW(:))-sum(tempRoi(:)))/sum(tempRoi(:))>.1
        error('problem with ROI construction')
    end
    
    
    % save sub-roi mask
    imwrite(tempRoi,[currentRoiAnDir filesep 'roiMask.tif']);
    save([currentRoiAnDir filesep 'roiYX'],'roiYX');
    
    % extract detected features falling in current ROI from the full set
    fNames=fieldnames(movieInfoOld);
    for iFrame=1:size(movieInfoOld,1)
        % look for feature coordinates within the polygon
        [inIdx,onIdx]=inpolygon(movieInfoOld(iFrame,1).xCoord(:,1),movieInfoOld(iFrame,1).yCoord(:,1),roiYX(:,2),roiYX(:,1));
        inPolyIdx=find(inIdx);
        % assign new movieInfo structure with only those features
        for iName=1:length(fNames)
            movieInfo(iFrame,1).(fNames{iName})=movieInfoOld(iFrame,1).(fNames{iName})(inPolyIdx,1);
        end    
    end
    % make feat directory for sub-roi and save movieInfo
    featDir=[currentRoiAnDir filesep 'feat']; mkdir(featDir);
    save([featDir filesep 'movieInfo'],'movieInfo')    

    
    % fill in new matrix of tracks
    tracksFinal=nan(size(trackedFeatureInfo));
    for iFrame=1:numTimePoints
        % find coordinate from tracks that are within polygon
        [inIdx,onIdx]=inpolygon(tracks.xCoord(:,iFrame),tracks.yCoord(:,iFrame),roiYX(:,2),roiYX(:,1));
        [inPolyIdx]=repmat(find(inIdx),[1 8]); % rows for points to add
        nPtsInside=size(inPolyIdx,1);

        col=repmat([sC(iFrame):eC(iFrame)],[nPtsInside,1]); % particular columns corresponding to frame
        
        tracksFinal(inPolyIdx(:),col(:))=trackedFeatureInfo(inPolyIdx(:),col(:));

    end
    % make track directory for sub-roi and save tracks matrix
    trackDir=[currentRoiAnDir filesep 'track']; mkdir(trackDir);
    save([trackDir filesep 'subTracks'],'tracksFinal')

    % ask user whether to not to select another sub-roi
    reply=input('Do you want to select another ROI? y/n [n]: ','s');
    if lower(reply)=='y'
        makeNewROI=1; % user said yes; make another one
        roiCount=roiCount+1; % counter for current condition rois
    else
        makeNewROI=0; % assume no; we're done
    end
end % while makeNewROI==1 && roiCount<10

% add a number to center of each sub-roi to show which region is which
imshow(img2show)
for iRoi=1:roiCount
    % make weighted mask using distance transform to find position where text should go 
    weightedRoi=bwdist(swapMaskValues(labelMatrix==iRoi));
    [r,c]=find(weightedRoi==max(weightedRoi(:)));      
    hold on
    text(c(1),r(1), num2str(iRoi),'color','r')
end

% save composite image and label matrix
frame = getframe(gca);
[I,map] = frame2im(frame);
imwrite(I,[pwd filesep 'sub-ROIs.tif'],'tif')
save('labelMatrix','labelMatrix');

close all

cd(homeDir)


function [img2show]=addMaskInColor(img,roiMask,c)
%subfunction to add new polygon outline to composite image
temp=double(bwmorph(roiMask,'remove'));
borderIdx=find(temp);
nPix=numel(roiMask);

img2show=img;
img2show(borderIdx)=c(1);
img2show(borderIdx+nPix)=c(2);
img2show(borderIdx+2*nPix)=c(3);

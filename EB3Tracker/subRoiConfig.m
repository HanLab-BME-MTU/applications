function subRoiConfig
% choose sub-regions of interest and get tracking data for them

% get project directory to work on
roiDir=uigetdir(pwd,'Please select ROI directory to use for making sub-ROIs');

% get images directory
homeDir=pwd;
cd(roiDir)
cd ..
imDir=[pwd filesep 'images'];


% load roiMask used to make ROI_x
roiMask=imread([roiDir filesep 'roiMask.tif']);


% load detected features
movieInfo=load([roiDir filesep 'feat' filesep 'movieInfo']);
movieInfoOld=movieInfo.movieInfo;

% load tracking data
% load tracksFinal (tracking result)
trackDir = [roiDir filesep 'track'];
[listOfFiles]=searchFiles('.mat',[],trackDir,0);
load([listOfFiles{1,2} filesep listOfFiles{1,1}])
% convert tracksFinal to matrix
[trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal);
clear trackedFeatureIndx

%get number of tracks and number of time points
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints=numTimePoints./8;
eC=8*[1:numTimePoints]';
sC=eC-7;

% save misc info for output

tracks.xCoord = trackedFeatureInfo(:,1:8:end); 
tracks.yCoord = trackedFeatureInfo(:,2:8:end);








% get list and number of images
[listOfImages]=searchFiles('.tif',[],imDir,0);

% read first image and make RGB
img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
img=(img-min(img(:)))./(max(img(:))-min(img(:)));
img=repmat(img,[1,1,3]);


% make sub-roi directory under ROI_x directory and go there
subRoiDir=[roiDir filesep 'subROIs'];
if isdir(subRoiDir)
    rmdir(subRoiDir,'s');
end
mkdir(subRoiDir);
cd(subRoiDir)


cMap=jet(10);


[img2show]=addMaskInColor(img,roiMask,cMap(10,:));

roiCount=1; % counter for rois for current project
makeNewROI=1; % flag for making new roi

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
    
    
    tempRoi=tempRoi & roiMask;
    roiMask=roiMask-tempRoi;
    labelMatrix(tempRoi)=roiCount;
   
    [img2show]=addMaskInColor(img2show,tempRoi,cMap(roiCount,:));
    
    roiYX=[polyYcoord polyXcoord; polyYcoord(1) polyXcoord(1)];

    % save sub-roi mask
    imwrite(tempRoi,[currentRoiAnDir filesep 'roiMask.tif']);
    
    
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

    
    % fill in new matrix
    tracksFinal=nan(size(trackedFeatureInfo));
    for iFrame=1:numTimePoints
        % find coordinate from tracks that are within polygon
        [inIdx,onIdx]=inpolygon(tracks.xCoord(:,iFrame),tracks.yCoord(:,iFrame),roiYX(:,2),roiYX(:,1));
        [inPolyIdx]=repmat(find(inIdx),[1 8]);
        nPtsInside=size(inPolyIdx,1);

        col=repmat([sC(iFrame):eC(iFrame)],[nPtsInside,1]);
        idx2add=repmat(numTracks.*[0:7],[nPtsInside 1]);
        allIdx=inPolyIdx+idx2add;
        tracksFinal(inPolyIdx(:),col(:))=trackedFeatureInfo(inPolyIdx(:),col(:));
    end
    trackDir=[currentRoiAnDir filesep 'track']; mkdir(trackDir);
    save([trackDir filesep 'subTracks'],'tracksFinal')
    
    reply=input('Do you want to select another ROI? y/n [n]: ','s');

    if lower(reply)=='y'
        makeNewROI=1; % user said yes; make another one
        roiCount=roiCount+1; % counter for current condition rois
    else
        makeNewROI=0; % assume no; we're done
    end
end % while makeNewROI==1 && roiCount<10

% add text to show which region is which
s = regionprops(labelMatrix, 'centroid');
centroids = cat(1, s.Centroid);
imshow(img2show)
for iRoi=1:roiCount
    hold on
    text(centroids(iRoi,1), centroids(iRoi,2), num2str(iRoi),'color','w')
end

frame = getframe(gca);
[I,map] = frame2im(frame);
imwrite(I,[pwd filesep 'sub-ROIs.tif'],'tif')
save('labelMatrix','labelMatrix');

close all

cd(homeDir)


function [img2show]=addMaskInColor(img,roiMask,c)

temp=double(bwmorph(roiMask,'remove'));
borderIdx=find(temp);
nPix=numel(roiMask);

img2show=img;
img2show(borderIdx)=c(1);
img2show(borderIdx+nPix)=c(2);
img2show(borderIdx+2*nPix)=c(3);

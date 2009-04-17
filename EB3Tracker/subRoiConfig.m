function subRoiConfig(projData)
% SUBROICONFIG allows user to choose sub-rois and get indices of subtracks
% which start and end in each sub-roi.

homeDir=pwd;

% load projData
if nargin<1 || isempty(projData)
    [fileName,pathName]=uigetfile('*.mat','Please select projData from META directory');
    if strcmp(fileName,0)
        return
    end
    projData=load([pathName filesep fileName]);
    projData=projData.projData;
end
anDir=formatPath(projData.anDir);
imDir=formatPath(projData.imDir);

% load roiYX and roiMask
roiYX=load([anDir filesep 'roiYX.mat']);
roiYX=roiYX.roiYX;
roiMask=imread([anDir filesep 'roiMask.tif']);

% get list of images, read first image, and make RGB (gray)
[listOfImages]=searchFiles('.tif',[],imDir,0);
img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
[imL,imW]=size(img(:,:,1));
img=(img-min(img(:)))./(max(img(:))-min(img(:)));
img=repmat(img,[1,1,3]);

% make sub-roi directory under roi_x directory and go there
subanDir=[anDir filesep 'subROIs'];
if isdir(subanDir)
    rmdir(subanDir,'s');
end
mkdir(subanDir);
cd(subanDir)

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

    % get coordinates of vertices of whole cell (ie all pixels of polygon boundary)
    [y1,x1]=ind2sub([imL,imW],find(roiMask,1)); % first pixel on boundary
    roiYXcell = bwtraceboundary(roiMask,[y1,x1],'N'); % get all pixels on boundary
    
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
    
    % find all subtracks starting within roi
    temp=projData.nTrack_start_end_velMicPerMin_class_lifetime;
    c=sub2ind(size(projData.xCoord),temp(:,1),temp(:,2));
    x=projData.xCoord(c);
    y=projData.yCoord(c);
    [inIdx,onIdx]=inpolygon(x,y,roiYX(:,2),roiYX(:,1));
    subtracksStartingIN=find(inIdx);
    save([currentRoiAnDir filesep 'subtracksStartingIN'],'subtracksStartingIN')
    
    % find all subtracks ending within roi
    c=sub2ind(size(projData.xCoord),temp(:,1),temp(:,3));
    x=projData.xCoord(c);
    y=projData.yCoord(c);
    [inIdx,onIdx]=inpolygon(x,y,roiYX(:,2),roiYX(:,1));
    subtracksEndingIN=find(inIdx);
    save([currentRoiAnDir filesep 'subtracksEndingIN'],'subtracksEndingIN')
    
    % ask user whether to not to select another sub-roi
    reply = questdlg('Do you want to select another ROI?');
    if strcmpi(reply,'yes')
        makeNewROI=1; % user said yes; make another one
        roiCount=roiCount+1; % counter for current condition rois
    else
        makeNewROI=0; % assume no; we're done
    end
end % while makeNewROI==1 && roiCount<10

% plot using vector graphics of boundaries and save as figure and tif
% add a number to center of each sub-roi to show which region is which
imshow(img)
hold on
% plot original roi outline
plot(roiYXcell(:,2),roiYXcell(:,1),'w');
for iRoi=1:roiCount
    % make weighted mask using distance transform to find position where text should go 
    weightedRoi=bwdist(swapMaskValues(labelMatrix==iRoi));
    [r,c]=find(weightedRoi==max(weightedRoi(:)));      
    text(c(1),r(1), num2str(iRoi),'color','r')
    % load sub-roi boundaries and plot outline
    currentRoiAnDir=[pwd filesep 'sub_' num2str(iRoi)];
    roiYX=load([currentRoiAnDir filesep 'roiYX']);
    roiYX=roiYX.roiYX;
    plot(roiYX(:,2),roiYX(:,1),'Color',cMap(iRoi,:));
end

% save composite image and label matrix
saveas(gcf,[pwd filesep 'sub-ROIs.fig'])
frame = getframe(gca);
[I,map] = frame2im(frame);
imwrite(I,[pwd filesep 'sub-ROIs.tif'],'tif')
save('labelMatrix','labelMatrix');

close all

cd(homeDir)


function [img2show]=addMaskInColor(img,roiMask,c)
%subfunction to add new polygon outline to composite image - this is needed
%because you can't pass vector graphics info to roipoly function, and we
%want to be able to visualize the regions that have already been selected.
temp=double(bwmorph(roiMask,'remove'));
borderIdx=find(temp);
nPix=numel(roiMask);

img2show=img;
img2show(borderIdx)=c(1);
img2show(borderIdx+nPix)=c(2);
img2show(borderIdx+2*nPix)=c(3);

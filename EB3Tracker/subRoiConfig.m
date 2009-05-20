function subRoiConfig(sourceProjData,fractionFromEdge)
% SUBROICONFIG allows user to choose sub-rois and get indices of subtracks
% which start and end in each sub-roi.

homeDir=pwd;
warningState=warning('query', 'all');

%warning on MATLAB:divideByZero

% load projData
if nargin<1 || isempty(sourceProjData)
    [fileName,roiMetaPathName]=uigetfile('*.mat','Please select projData from META directory');
    if strcmp(fileName,0)
        return
    end
    sourceProjData=load([roiMetaPathName filesep fileName]);
    sourceProjData=sourceProjData.projData;
end
anDir=formatPath(sourceProjData.anDir);
imDir=formatPath(sourceProjData.imDir);

if nargin<1 || isempty(fractionFromEdge)
    fractionFromEdge=[];
end

% load roiYX and roiMask
roiYX=load([anDir filesep 'roiYX.mat']);
roiYX=roiYX.roiYX;
roiMask=imread([anDir filesep 'roiMask.tif']);


roiArea=sum(roiMask(:));

% get list of images, read first image, and make RGB (gray)
[listOfImages]=searchFiles('.tif',[],imDir,0);
img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
[imL,imW]=size(img(:,:,1));
img1=(img-min(img(:)))./(max(img(:))-min(img(:)));
img=repmat(img1,[1,1,3]);

if ~isempty(fractionFromEdge)

    % get roi centroid
    stats=regionprops(bwlabel(roiMask),'centroid');
    centerRoiYX=stats.Centroid(2:-1:1);

    % make distance transform
    weightedRoi=bwdist(swapMaskValues(roiMask));
    % innerMask's periphery is fractionFromEdge away from the original roi boundary
    innerMask=weightedRoi>fractionFromEdge.*max(weightedRoi(:));

    % get list of all the pixels on innerRoi's boundary
    [y1,x1]=ind2sub([imL,imW],find(innerMask,1));
    if fractionFromEdge<1
        innerMaskYX = bwtraceboundary(innerMask,[y1,x1],'N');
    else
        innerMaskYX=[nan nan];
    end

    % get list of all the pixels on innerRoi's boundary
    [y1,x1]=ind2sub([imL,imW],find(roiMask,1));
    roiMaskYX = bwtraceboundary(roiMask,[y1,x1],'N');

    figure
    imshow(roiMask.*img1,[])
    hold on
    axis equal
    plot(innerMaskYX(:,2),innerMaskYX(:,1))
    plot(roiMaskYX(:,2),roiMaskYX(:,1))
    h=msgbox('Draw a line and double-click when finished','help');
    uiwait(h);
    h=imline;
    position = wait(h);
    close(gcf)

    % position of the ends of the user-chosen line
    lineEndsYX=position(:,2:-1:1);

    [xAll,yAll]=meshgrid(1:imW,1:imL);
    if abs(lineEndsYX(1,1)-lineEndsYX(2,1))<1 % y's are the same, user chose horizontal line
        r11=(yAll<=lineEndsYX(1,1));
        r12=(yAll> lineEndsYX(1,1));
        r21=(xAll<=centerRoiYX(1,2));
        r22=(xAll> centerRoiYX(1,2));
    elseif abs(lineEndsYX(1,2)-lineEndsYX(2,2))<1 % x's are the same, user chose vertical line
        r11=(xAll<=lineEndsYX(1,2));
        r12=(xAll> lineEndsYX(1,2));
        r21=(yAll<=centerRoiYX(1,1));
        r22=(yAll> centerRoiYX(1,1));
    else
        % user-chosen line slope and y-intercept
        m1=diff(lineEndsYX(:,1))/diff(lineEndsYX(:,2));
        b1=lineEndsYX(1,1)-m1*lineEndsYX(1,2);
        % perpendicular line going through roi's centroid
        m2=-1/m1;
        b2=centerRoiYX(1,1)-m2*centerRoiYX(1,2);
        % y-coordinates of both lines across all the x-pixels
        yLine1=repmat(m1.*[1:imW]+b1,[imL,1]);
        yLine2=repmat(m2.*[1:imW]+b2,[imL,1]);
        % divide the image into two parts on either side of line 1
        r11=(yAll<=yLine1);
        r12=(yAll> yLine1);
        % divide the image into two parts on either side of line 2
        r21=(yAll<=yLine2);
        r22=(yAll> yLine2);

    end
    
    if fractionFromEdge<1
        roiSet=zeros(imL,imW,5);
        roiSet(:,:,1)=innerMask;
        roiSet(:,:,2)=r11 & r21 & (roiMask-innerMask);
        roiSet(:,:,3)=r11 & r22 & (roiMask-innerMask);
        roiSet(:,:,4)=r12 & r21 & (roiMask-innerMask);
        roiSet(:,:,5)=r12 & r22 & (roiMask-innerMask);
    else
        roiSet=zeros(imL,imW,4);
        roiSet(:,:,1)=r11 & r21 & (roiMask-innerMask);
        roiSet(:,:,2)=r11 & r22 & (roiMask-innerMask);
        roiSet(:,:,3)=r12 & r21 & (roiMask-innerMask);
        roiSet(:,:,4)=r12 & r22 & (roiMask-innerMask);
    end

end




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
    mkdir([currentRoiAnDir filesep 'meta']);
    
    if ~isempty(fractionFromEdge)
        tempRoi=roiSet(:,:,roiCount);
    else
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
    end

    % get coordinates of vertices of whole cell (ie all pixels of polygon boundary)
    [y1,x1]=ind2sub([imL,imW],find(roiMask,1)); % first pixel on boundary
    roiYXcell = bwtraceboundary(roiMask,[y1,x1],'N'); % get all pixels on boundary

    % get intersection with max region in the cell
    tempRoi=tempRoi & roiMask;
    % shrink max region in the cell for next round by excluding current roi
    roiMask=roiMask-tempRoi;
    % fill in label matrix
    labelMatrix(tempRoi)=roiCount;

    percentRoiArea=sum(tempRoi(:))/roiArea;
    
    allArea(roiCount)=percentRoiArea;

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
    save([currentRoiAnDir filesep 'percentRoiArea'],'percentRoiArea');

    % assign big matrix with shorter name
    aT=sourceProjData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;

    % find all subtracks where the midpoint is within the subroi
    c=sub2ind(size(sourceProjData.xCoord),aT(:,1),round(mean([aT(:,2),aT(:,3)],2)));
    x=sourceProjData.xCoord(c);
    y=sourceProjData.yCoord(c);
    [inIdx,onIdx]=inpolygon(x,y,roiYX(:,2),roiYX(:,1));
    subIdx=find(inIdx);
    save([currentRoiAnDir filesep 'subIdx'],'subIdx')

    % projData will have same format as sourceProjData but with only data for
    % relevant tracks
    projData=sourceProjData;

    % project name is now the sub-roi directory
    projData.anDir=currentRoiAnDir;

    % initialize matrices with all nans
    projData.xCoord=nan.*sourceProjData.xCoord;
    projData.yCoord=nan.*sourceProjData.yCoord;
    projData.featArea=nan.*sourceProjData.featArea;
    projData.featInt=nan.*sourceProjData.featInt;
    projData.frame2frameVel_micPerMin=nan.*sourceProjData.frame2frameVel_micPerMin;
    projData.segGapAvgVel_micPerMin=nan.*sourceProjData.segGapAvgVel_micPerMin;

    % reassign appropriate sections of the data for relevant subtracks
    for iSub=1:length(subIdx)
        row =aT(subIdx(iSub),1); % full track number
        cols=aT(subIdx(iSub),2):aT(subIdx(iSub),3); % relevant frames
        projData.xCoord(row,cols)   = sourceProjData.xCoord(row,cols);
        projData.yCoord(row,cols)   = sourceProjData.yCoord(row,cols);
        projData.featArea(row,cols) = sourceProjData.featArea(row,cols);
        projData.featInt(row,cols)  = sourceProjData.featInt(row,cols);
        projData.frame2frameVel_micPerMin(row,cols(1:end-1))= sourceProjData.frame2frameVel_micPerMin(row,cols(1:end-1));
        projData.segGapAvgVel_micPerMin(row,cols(1:end-1))  = sourceProjData.segGapAvgVel_micPerMin(row,cols(1:end-1));
    end

    % reassign matrix with only the subset of tracks in the subroi
    aT=aT(subIdx,:); % just use a shorter name to make it easier
    projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix=aT;

    % count the tracks that appear in the sub-roi
    % projData.numTracks=length(unique(aT(:,1)));

    % don't worry about these parameters for now...
    projData.frame2frameDispPix=[];
    projData.pair2pairDiffPix=[];
    projData.medNNdistWithinFramePix=[];
    projData.meanDisp2medianNNDistRatio=[];


    % figure out which track numbers contain a pause, catastrophe, or
    % unclassified gap
    projData.tracksWithPause       = unique(aT(aT(:,5)==2,1));
    projData.tracksWithCatastrophe = unique(aT(aT(:,5)==3,1));
    projData.tracksWithUnclassifed = unique(aT(aT(:,5)==4,1));

    % mean/std for growth speed (microns per minute)
    projData.typeStats.type1_mean_micPerMin = mean(aT(aT(:,5)==1,4));
    projData.typeStats.type1_std_micPerMin  =  std(aT(aT(:,5)==1,4));

    % probability of pausing is 1 over the average total time (in seconds) spent
    % growing prior to pause event
    pauseIdx=find(aT(:,5)==2);
    beforePauseIdx=pauseIdx-1;
    % if the first row is a pause, the growth phase before isn't in the
    % sub-roi.  get rid of this one.
    if beforePauseIdx(1)==0
        pauseIdx(1)=[];
        beforePauseIdx(1)=[];
    end
    % get rid of those indices corresponding to the growth phase of a
    % different track OR where the growth phase begins in the first frame
    remIdx=find(aT(pauseIdx,1)~=aT(beforePauseIdx,1) | aT(beforePauseIdx,2)==1);
    beforePauseIdx(remIdx)=[];

    if isempty(beforePauseIdx)
        projData.typeStats.type1_Ppause=NaN;
    else
        projData.typeStats.type1_Ppause=mean(1./(aT(beforePauseIdx,6).*sourceProjData.secPerFrame));
    end

    % probability of shrinking is 1 over the average total time (in seconds) spent
    % growing prior to shrinkage event
    shrinkIdx=find(aT(:,5)==3);
    beforeShrinkIdx=shrinkIdx-1;
    % if the first row is a catastrophe, the growth phase before isn't in the
    % sub-roi.  get rid of this one.
    if beforeShrinkIdx(1)==0
        shrinkIdx(1)=[];
        beforeShrinkIdx(1)=[];
    end

    % get rid of those indices corresponding to the growth phase of a
    % different track OR where the growth phase begins in the first frame
    remIdx=find(aT(shrinkIdx,1)~=aT(beforeShrinkIdx,1) | aT(beforeShrinkIdx,2)==1);
    beforeShrinkIdx(remIdx)=[];

    if isempty(beforeShrinkIdx)
        projData.typeStats.type1_Pcat=NaN;
    else
        projData.typeStats.type1_Pcat=mean(1./(aT(beforeShrinkIdx,6).*sourceProjData.secPerFrame));
    end

    % mean/std for pause speed (microns per minute)
    if ~isempty(aT(aT(:,5)==2,4))
        projData.typeStats.type2_mean_micPerMin = mean(aT(aT(:,5)==2,4));
        projData.typeStats.type2_std_micPerMin  =  std(aT(aT(:,5)==2,4));
    else
        projData.typeStats.type2_mean_micPerMin = NaN;
        projData.typeStats.type2_std_micPerMin  = NaN;
    end

    % mean/std for shrinkage speed (microns per minute)
    if ~isempty(aT(aT(:,5)==3,4))
        projData.typeStats.type3_mean_micPerMin = mean(aT(aT(:,5)==3,4));
        projData.typeStats.type3_std_micPerMin  =  std(aT(aT(:,5)==3,4));
    else
        projData.typeStats.type3_mean_micPerMin = NaN;
        projData.typeStats.type3_std_micPerMin  = NaN;
    end

    % mean/std for unclassified speed (microns per minute)
    if ~isempty(aT(aT(:,5)==4,4))
        projData.typeStats.type4_mean_micPerMin = mean(aT(aT(:,5)==4,4));
        projData.typeStats.type4_std_micPerMin  =  std(aT(aT(:,5)==4,4));
    else
        projData.typeStats.type4_mean_micPerMin = NaN;
        projData.typeStats.type4_std_micPerMin  = NaN;
    end


    save([currentRoiAnDir filesep 'meta' filesep 'projData'],'projData')

    if ~isempty(fractionFromEdge)
        if roiCount<size(roiSet,3)
            reply='yes';
        else
            reply='no';
        end
    else
        reply = questdlg('Do you want to select another ROI?');
    end
    if strcmpi(reply,'yes')
            makeNewROI=1; % user said yes; make another one
            roiCount=roiCount+1; % counter for current condition rois
        else
            makeNewROI=0; % assume no; we're done
        end
    end % while makeNewROI==1 && roiCount<10

% plot using vector graphics of boundaries and save as figure and tif
% add a number to center of each sub-roi to show which region is which
figure
imshow(img)
hold on
% plot original roi outline
plot(roiYXcell(:,2),roiYXcell(:,1),'w');
for iRoi=1:roiCount
    % make weighted mask using distance transform to find position where text should go
    weightedRoi=bwdist(swapMaskValues(labelMatrix==iRoi));
    [r,c]=find(weightedRoi==max(weightedRoi(:)));
    text(c(1),r(1), {['Sub-roi: ' num2str(iRoi)],['Fraction of area: ' sprintf('%3.2f',allArea(iRoi))]},'color','r');
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

cd(anDir)
cd ..
getProj(pwd);


cd(homeDir)
warning(warningState);


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

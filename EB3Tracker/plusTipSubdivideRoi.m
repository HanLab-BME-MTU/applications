function plusTipSubdivideRoi(sourceProjData,fractionFromEdge,savedROI)
% plusTipSubdivideRoi allows user to choose sub-regions of interest
%
% INPUT : sourceProjData        : projData file from any project
%         fractionFromEdge (OPT): decimal value representing position of a
%                                 central polygon edge from ROI-boundary to
%                                 the cell's center of mass.  Use 1 to
%                                 avoid the central polygon and only have
%                                 four quadrants.  Value should be in range
%                                 0-1. Upon running the function, you will
%                                 be asked to draw a line across the cell.
%                                 This provides one of the two axes used to
%                                 divide the cell into quadrants. If this
%                                 parameter is not included, the user can
%                                 pick his/her own regions.
%         savedROI              : either a BW mask or roiYX coordinates. if
%                                 empty, user picks a new ROI
%
% OUTPUT: sub-roi directories, a tiff image showing regional selections,
%         and a text file giving the speed/lifetime/displacement
%         distributions for growth.  a growth track is included in a
%         particular region if more than half the frames of that track's
%         existence fall within the region. merged tracks are used in the
%         case of reclassified fgaps.

doPlot=0;

homeDir=pwd;
warningState = warning;
warning('off','Images:initSize:adjustingMag')

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

if ~isempty(strfind(anDir,'sub'))
    msgbox('Cannot choose sub-ROIs from a sub-ROI')
    return
end


cd(anDir)
imDir=formatPath(sourceProjData.imDir);

if nargin<2 || isempty(fractionFromEdge)
    fractionFromEdge=[];
elseif fractionFromEdge<0 || fractionFromEdge>1
    msgbox('Fraction must be in 0-1 range')
    return   
end

% make sub-roi directory under roi_x directory and go there
subanDir=[anDir filesep 'subROIs'];
if ~isdir(subanDir)
    mkdir(subanDir);
end
cd(subanDir)

flag=0;
roiCount=1; % counter for rois for current project
while flag==0
    currentRoiAnDir=[pwd filesep 'sub_' num2str(roiCount)];
    if isdir(currentRoiAnDir)
        roiCount=roiCount+1;
    else
        flag=1;
    end
end
roiStart=roiCount;

% input to varycolor is black (we don't use it since the image is dark;
% make it one more than the number of ROIs (max is 9)
cMap=varycolor(10);


% get list of images, read first image, and make RGB (gray)
[listOfImages]=searchFiles('.tif',[],imDir,0);
img=double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]));
[imL,imW]=size(img(:,:,1));
img1=(img-min(img(:)))./(max(img(:))-min(img(:)));
img=repmat(img1,[1,1,3]);

% get figure size for parallel coordinates plot
% scrsz = get(0,'ScreenSize');
% screenW=scrsz(3);
% screenL=scrsz(4);
% magCoef=inf;
% xRange=imW;
% yRange=imL;
% maxMagCoefW = (0.8*screenW)/xRange;
% maxMagCoefL = (0.8*screenL)/yRange;
% if magCoef > min([maxMagCoefW; maxMagCoefL])
%     calcMagCoef = min([magCoef; maxMagCoefW; maxMagCoefL]);
% else
%     calcMagCoef = magCoef;
% end
% figW = (calcMagCoef*xRange);
% figL = (calcMagCoef*yRange);
% figPos=[round(screenW*(1-figW/screenW)/2) round(screenL*(1-figL/screenL)/2) figW figL];



% get roiYX and roiMask
if nargin<3 || isempty(savedROI)
    h=msgbox('First draw the full region of interest ROI','help');
    uiwait(h);
    roiMask=[];
    while isempty(roiMask)
        try
            figure %('Position',figPos)
            roiMask=roipoly(img);
            close

        catch
            h=msgbox('Please try again.','help');
            uiwait(h);
        end
    end
else
    if ~isequal(size(savedROI),size(img(:,:,1)))
        roiMask=roipoly(img,savedROI(:,2),savedROI(:,1));
    else
        roiMask=savedROI;
    end
end
roiArea=sum(roiMask(:));
% save roiMask and coordinates
imwrite(roiMask,[subanDir filesep 'fullRoiMask_' num2str(roiStart) '.tif']);

% string for number of files
strg1 = sprintf('%%.%dd',2);


% if the input includes the fraction from the cell edge parameter, then
% divide the cell into 5 sub-rois consisting of a central polygon and four
% quadrants around the periphery
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

    figure %('Position',figPos)
    imshow(roiMask.*img1,[])
    hold on
    axis equal
    plot(innerMaskYX(:,2),innerMaskYX(:,1))
    plot(roiMaskYX(:,2),roiMaskYX(:,1))
    h=msgbox('Draw a line across the cell and double-click when finished','help');
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

% set cell boundary to white in composite image
[img2show]=addMaskInColor(img,roiMask,[1 1 1]);

% store 1 in every pixel of sub_1, 2 in every pixel of sub_2, etc.
labelMatrix=zeros(size(roiMask));

% iterate til the user is finished
makeNewROI=1; % flag for making new roi
while makeNewROI==1
    % make new roi_n image/analysis directories
    %indxStr1 = sprintf(strg1,roiCount);
    indxStr1 = num2str(roiCount);
    currentRoiAnDir=[pwd filesep 'sub_' indxStr1];
    mkdir(currentRoiAnDir);
    mkdir([currentRoiAnDir filesep 'meta']);

    sourceFeatDir=[sourceProjData.anDir filesep 'feat'];
    mkdir([currentRoiAnDir filesep 'feat']);
    copyfile([anDir filesep 'feat' filesep 'movieInfo.mat'],[currentRoiAnDir filesep 'feat' filesep 'movieInfo.mat']);

    if ~isempty(fractionFromEdge)
        tempRoi=roiSet(:,:,roiCount-roiStart+1);
    else
        tempRoi=[];
        while isempty(tempRoi)
            try
                if makeNewROI==1
                    h=msgbox('Draw a sub-ROI','help');
                    uiwait(h);
                end
                figure %('Position',figPos)
                [tempRoi,polyXcoord,polyYcoord]=roipoly(img2show);
            catch
                disp('Please try again.')
            end
        end
        close(gcf)
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

    totalAreaPixels=sum(tempRoi(:));
    percentRoiArea=100*(totalAreaPixels/roiArea);

    allArea(roiCount,1:2)=[totalAreaPixels percentRoiArea];

    % add the current roi to the composite image
    [img2show]=addMaskInColor(img2show,tempRoi,cMap(max(1,mod(roiCount,10)),:));

    % get coordinates of vertices (ie all pixels of polygon boundary)
    [y1,x1]=ind2sub([imL,imW],find(tempRoi,1)); % first pixel on boundary
    roiYX = bwtraceboundary(tempRoi,[y1,x1],'N'); % get all pixels on boundary
    % test to make sure roi can be reproduced ok - assume that new polygon
    % doesn't differ in area more than 10% of tempRoi
    [resultBW]=roipoly(imL,imW,roiYX(:,2),roiYX(:,1));
    if abs(sum(resultBW(:))-sum(tempRoi(:)))/sum(tempRoi(:))>.1
        error('problem with ROI construction')
    end

    % projData will have same format as sourceProjData but with only data for
    % relevant tracks
    projData=sourceProjData;
    projData.anDir=currentRoiAnDir;

    [aTmerge,aTreclass,dataMatCrpSecMic]=plusTipMergeSubtracks(projData);

    % only retain growth subtracks from aTmerge
    aTmerge(aTmerge(:,5)~=1,:)=[];
    % only retain subtracks where the midpoint is within the subroi
    c=sub2ind(size(sourceProjData.xCoord),aTmerge(:,1),round(mean([aTmerge(:,2),aTmerge(:,3)],2)));
    x=sourceProjData.xCoord(c);
    y=sourceProjData.yCoord(c);
    [inIdx,onIdx]=inpolygon(x,y,roiYX(:,2),roiYX(:,1));
    subIdx=find(inIdx);
    aTmerge=aTmerge(subIdx,:);

    % only retain growth subtracks from dataMatCrpSecMic
    dataMatCrpSecMic(dataMatCrpSecMic(:,5)~=1,:)=[];
    % only retain subtracks where the midpoint is within the subroi
    c=sub2ind(size(sourceProjData.xCoord),dataMatCrpSecMic(:,1),round(mean([dataMatCrpSecMic(:,2),dataMatCrpSecMic(:,3)],2)));
    x=sourceProjData.xCoord(c);
    y=sourceProjData.yCoord(c);
    [inIdx,onIdx]=inpolygon(x,y,roiYX(:,2),roiYX(:,1));
    subIdx=find(inIdx);
    dataMatCrpSecMic=dataMatCrpSecMic(subIdx,:);


    % new number of tracks
    projData.numTracks=length(unique(aTmerge(:,1)));

    % keep only the coordinates, speeds, etc. corresponding to tracks remaining
    projData.xCoord=nan(size(sourceProjData.xCoord));
    projData.yCoord=nan(size(sourceProjData.yCoord));
    projData.featArea=nan(size(sourceProjData.featArea));
    projData.featInt=nan(size(sourceProjData.featInt));
    projData.frame2frameVel_micPerMin=nan(size(sourceProjData.frame2frameVel_micPerMin));
    projData.segGapAvgVel_micPerMin=nan(size(sourceProjData.segGapAvgVel_micPerMin));
    for iSub=1:size(aTmerge,1)
        projData.xCoord(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3))=sourceProjData.xCoord(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3));
        projData.yCoord(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3))=sourceProjData.yCoord(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3));

        projData.featArea(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3))=sourceProjData.featArea(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3));
        projData.featInt(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3))=sourceProjData.featInt(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3));

        projData.frame2frameVel_micPerMin(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3)-1)=sourceProjData.frame2frameVel_micPerMin(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3)-1);
        projData.segGapAvgVel_micPerMin(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3)-1)=sourceProjData.segGapAvgVel_micPerMin(aTmerge(iSub,1),aTmerge(iSub,2):aTmerge(iSub,3)-1);
    end


    % get frame-to-frame displacement for growth only (not forward/backward gaps)
    frame2frameDispPix=sqrt(diff(projData.xCoord,1,2).^2+diff(projData.yCoord,1,2).^2);
    % get rid of NaNs and linearize the vector
    projData.frame2frameDispPix=frame2frameDispPix(~isnan(frame2frameDispPix(:)));

    % get change in velocity between frame *pairs* for segments only
    pair2pairDiffPix=diff(frame2frameDispPix,1,2);
    % get rid of NaNs and linearize the vector
    projData.pair2pairDiffPix=pair2pairDiffPix(~isnan(pair2pairDiffPix(:)));
    % std (microns/min) of delta growthSpeed btw frames
    projData.pair2pairDiffMicPerMinStd=std(pixPerFrame2umPerMin(projData.pair2pairDiffPix,projData.secPerFrame,projData.pixSizeNm));

    projData.medNNdistWithinFramePix=NaN;
    projData.meanDisp2medianNNDistRatio=NaN;

    % there are no track numbers that contain an fgap or bgap
    projData.percentFgapsReclass=NaN;
    projData.percentBgapsReclass=NaN;
    projData.tracksWithFgap = NaN;
    projData.tracksWithBgap = NaN;

    % calculate stats using the matrix where beginning/end data has been
    % removed. M records speeds (microns/min), lifetimes (sec), and
    % displacements (microns) for growths, fgaps,and bgaps.
    [projData.stats,M]=plusTipDynamParam(dataMatCrpSecMic);

    projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix=aTmerge;

    % save sub-roi mask and coordinates
    imwrite(tempRoi,[currentRoiAnDir filesep 'roiMask.tif']);
    save([currentRoiAnDir filesep 'roiYX'],'roiYX');
    % save each projData in its own directory
    save([currentRoiAnDir filesep 'meta' filesep 'projData'],'projData')
    save([currentRoiAnDir filesep 'subRoiInfo'],'totalAreaPixels','percentRoiArea');

    % write out speed/lifetime/displacement distributions into a text file
    dlmwrite([currentRoiAnDir filesep 'gs_fs_bs_gl_fl_bl_gd_fd_bd.txt'], M, 'precision', 3,'delimiter', '\t','newline', 'pc');

    if doPlot==1

        figure;
        hist(sourceProjData.frame2frameDispPix,50)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r')
        hold on;
        hist(projData.frame2frameDispPix,50)
        title('growth phase frame-to-frame displacement (pixels)')
        legend('ROI','Sub-ROI')

        figure;
        hist(sourceProjData.pair2pairDiffPix,50)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r')
        hold on;
        hist(projData.pair2pairDiffPix,50)
        title('growth phase frame-pair-to-frame-pair difference (pixels)')
        legend('ROI','Sub-ROI')

    end

    if ~isempty(fractionFromEdge)
        if roiCount<size(roiSet,3)+roiStart-1
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
figure %('Position',figPos)
imshow(img)
hold on
% plot original roi outline
plot(roiYXcell(:,2),roiYXcell(:,1),'w');
for iRoi=roiStart:roiCount
    % make weighted mask using distance transform to find position where text should go
    weightedRoi=bwdist(swapMaskValues(labelMatrix==iRoi));
    [r,c]=find(weightedRoi==max(weightedRoi(:)));
    %indxStr1 = sprintf(strg1,iRoi);
    indxStr1 = num2str(iRoi);
    text(c(1),r(1), {['Sub-ROI: ' indxStr1],['Total pixels: ' sprintf('%3.2f',allArea(iRoi,1))],[sprintf('%3.2f',allArea(iRoi,2)) ' %']},'color','r');
    % load sub-roi boundaries and plot outline
    currentRoiAnDir=[pwd filesep 'sub_' indxStr1];
    roiYX=load([currentRoiAnDir filesep 'roiYX']);
    roiYX=roiYX.roiYX;
    plot(roiYX(:,2),roiYX(:,1),'Color',cMap(max(1,mod(iRoi,10)),:));
end

% save composite image and label matrix

frame = getframe(gca);
[I,map] = frame2im(frame);
% indxStr1 = sprintf(strg1,roiStart);
% indxStr2 = sprintf(strg1,roiCount);
indxStr1 = num2str(roiStart);
indxStr2 = num2str(roiCount);
imwrite(I,[pwd filesep 'subROIs_' indxStr1 '_' indxStr2 '.tif'])
saveas(gcf,[pwd filesep 'subROIs_' indxStr1 '_' indxStr2 '.fig'])
save(['labelMatrix_' indxStr1 '_' indxStr2],'labelMatrix');

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

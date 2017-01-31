function h2 = showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,IDtoInspect, gPath,additionalName)
% h2 = showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,IDtoInspect, gPath,additionalName)
% This function shows big picture, montage, and time series of fluorescence
% intensity and traction.
% Sangyoon Han, Jun 2015
imgWidth = size(imgMap,2);
imgHeight = size(imgMap,1);
curEndFrame=curTrack.endingFrameExtra;
curStartFrame = curTrack.startingFrameExtra;
curEndFrameEE=curTrack.endingFrameExtraExtra;
curStartFrameEE = curTrack.startingFrameExtraExtra;
curFrameRange= curStartFrameEE:curEndFrameEE;
chosenStartFrame = curStartFrameEE;
chosenEndFrame = curEndFrameEE;

pixSize = MD.pixelSize_; % nm/pixel
tInterval = MD.timeInterval_; % time interval in sec
scaleBar = 1; %micron

% curFrameRangeEE= curStartFrameEE:curEndFrameEE;
chosenFRange = curFrameRange;
[~,peakFrame] = max(curTrack.ampTotal(chosenFRange));
peakFrame = chosenFRange(peakFrame);
if ~isempty(curTrack.intenPeakness) && curTrack.intenPeakness
    peakFrame = curTrack.intenPeakFrame;
end

% ROI
%     r_pix = 5; % half the width of the ROI
% this r_pix was shown to be too small if the maturing adhesion slides.
% r_pix should include all of the track trace ...
maxX = nanmax(curTrack.xCoord)-nanmin(curTrack.xCoord);
maxY =  nanmax(curTrack.yCoord)-nanmin(curTrack.yCoord);
r_pix = ceil(max(max(maxX,maxY)/2+7,10));
meanX = round(nanmean(curTrack.xCoord));
meanY = round(nanmean(curTrack.yCoord));
bLeft = max(1,meanX-r_pix);
bRight = min(imgWidth,meanX+r_pix);
bBottom = max(1,meanY-r_pix);
bTop = min(imgHeight,meanY+r_pix);
% ROI for big images
r_pixBig = 90; % half the width of the ROI
bLeftBig = max(1,meanX-r_pixBig);
bRightBig = min(imgWidth,meanX+r_pixBig);
bBottomBig = max(1,meanY-r_pixBig);
bTopBig = min(imgHeight,meanY+r_pixBig);
% make a new figure
h2=figure;
figWidth = 600;
figHeight = 675;
%     set(h2,'Position',[100,10,figWidth,figHeight]),title(['ID:' num2str(IDtoInspect)])
set(h2,'Position',[100,200,figWidth,figHeight])
set(h2,'PaperPositionMode','auto')
set(h2,'Resize','off')
set(h2,'Renderer','painters')
set(h2,'InvertHardcopy','off')

genFontSize = 7;
genFontUnit = 'points';
% Overall shots 
imgBigXLim = [bLeftBig, bRightBig];
imgBigYLim = [bBottomBig, bTopBig];
marginX = 5/figWidth; %normalized margin width
marginY = 5/figHeight; %normalized margin width
mgScale = 8;
% ax1: initial point
% it would be better to show frame that actually shows adhesion, so we
% will use frame that is mid-point between initial and peak
midStartFrame = ceil((chosenStartFrame+peakFrame)/2);
ax1=axes('Position',[marginX, 490/figHeight, 175/figWidth,175/figHeight]);
%get the dynamic range
cropImg = imgMap(bBottom:bTop,bLeft:bRight,chosenStartFrame:chosenEndFrame);
cropImgSmaller = imgMap(bBottom+7:bTop-7,bLeft+7:bRight-7,chosenStartFrame:chosenEndFrame);
maxInt = max(cropImgSmaller(:));
minInt = min(cropImgSmaller(:));
imshow(imgMap(:,:,midStartFrame),[minInt maxInt]), hold on
% have to show a box indicating a region of interest
rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
set(ax1,'XLim',imgBigXLim,'YLim',imgBigYLim)
line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
    [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
midStartFrameAdjusted = midStartFrame - chosenStartFrame; % in frames
txtAx1= text(bRightBig-40, bBottomBig+12,[num2str(midStartFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
set(txtAx1,'horizontalAlignment','right')
set(txtAx1,'position',[bRightBig-mgScale,bBottomBig+12])
set(findobj(ax1,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)
% ax2: peak point
ax2=subplot('Position',[2*marginX+175/figWidth, 490/figHeight, 175/figWidth, 175/figHeight]);
imshow(imgMap(:,:,peakFrame),[minInt maxInt]), hold on
rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
set(ax2,'XLim',imgBigXLim,'YLim',imgBigYLim)
line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
    [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
peakFrameAdjusted = peakFrame - chosenStartFrame; % in frames
txtAx2= text(bRightBig-40, bBottomBig+12,[num2str(peakFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
set(txtAx2,'horizontalAlignment','right')
set(txtAx2,'position',[bRightBig-mgScale,bBottomBig+12])
set(findobj(ax2,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)
% ax3: end point
% it would be better to show frame that actually shows adhesion, so we
% will use frame that is mid-point between initial and peak
midEndFrame = ceil((chosenEndFrame+peakFrame)/2);
ax3=subplot('Position',[3*marginX+2*175/figWidth, 490/figHeight, 175/figWidth, 175/figHeight]);
imshow(imgMap(:,:,midEndFrame),[minInt maxInt]), hold on
rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
set(ax3,'XLim',imgBigXLim,'YLim',imgBigYLim)
line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
    [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
midEndFrameAdjusted = midEndFrame - chosenStartFrame; % in frames
txtAx3= text(bRightBig-40, bBottomBig+12,[num2str(midEndFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
set(txtAx3,'horizontalAlignment','right')
set(txtAx3,'position',[bRightBig-mgScale,bBottomBig+12])
set(findobj(ax3,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)
% ax4: initial point with traction
%     tMap2 = tMap(:,:,[chosenStartFrame,peakFrame,chosenEndFrame]);
tMap2 = tMap(:,:,[midStartFrame,peakFrame,midEndFrame]);
tMapROI = tMap2(bBottomBig:bTopBig,bLeftBig:bRightBig,:);
tCropMax = max(curTrack.forceMag(chosenStartFrame:chosenEndFrame))*1.1;
tCropMin = min(curTrack.forceMag(chosenStartFrame:chosenEndFrame));
tmax= min(tCropMax*3,max(tMapROI(:))*0.9);
tmin= tCropMin*0.1;
% tmax= max(tMapROI(:))*0.9;
% tmin= min(tMapROI(:));

ax4=subplot('Position',[marginX, 250/figHeight+marginY, 175/figWidth,175/figHeight]);
imshow(tMap2(:,:,1),[tmin tmax]), hold on
rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
set(ax4,'XLim',imgBigXLim,'YLim',imgBigYLim)
line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
    [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
txtScaleAx4 = text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
set(txtScaleAx4,'FontUnits',genFontUnit)
set(txtScaleAx4,'FontSize',genFontSize)
txtAx4= text(bRightBig-40, bBottomBig+12,[num2str(midStartFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
set(txtAx4,'horizontalAlignment','right')
set(txtAx4,'position',[bRightBig-mgScale,bBottomBig+12])

set(findobj(ax4,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)
% ax5: fluorescence peak point with traction
ax5=subplot('Position',[2*marginX+175/figWidth, 250/figHeight+marginY, 175/figWidth,175/figHeight]);
imshow(tMap2(:,:,2),[tmin tmax]), hold on
rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
set(ax5,'XLim',imgBigXLim,'YLim',imgBigYLim)
line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
    [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
txtAx5= text(bRightBig-40, bBottomBig+12,[num2str(peakFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
set(txtAx5,'horizontalAlignment','right')
set(txtAx5,'position',[bRightBig-mgScale,bBottomBig+12])
set(findobj(ax5,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)
% ax6: end point with traction
ax6=subplot('Position',[3*marginX+2*175/figWidth, 250/figHeight+marginY, 175/figWidth,175/figHeight]);
imshow(tMap2(:,:,3),[tmin tmax]), hold on
rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
set(ax6,'XLim',imgBigXLim,'YLim',imgBigYLim)
line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
    [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
txtAx6= text(bRightBig-40, bBottomBig+12,[num2str(midEndFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
set(txtAx6,'horizontalAlignment','right')
set(txtAx6,'position',[bRightBig-mgScale,bBottomBig+12])
set(findobj(ax6,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)
%Montage
% cropImg = imgMap(bBottom:bTop,bLeft:bRight,chosenStartFrame:chosenEndFrame);
cropTmap = tMap(bBottom:bTop,bLeft:bRight,chosenStartFrame:chosenEndFrame);
% make 4D array to use montage
cropImg = reshape(cropImg,size(cropImg,1),size(cropImg,2),1,size(cropImg,3));
cropTmap = reshape(cropTmap,size(cropTmap,1),size(cropTmap,2),1,size(cropTmap,3));
axes('Position',[marginX, 430/figHeight+marginY, 540/figWidth-marginX,50/figHeight]);
% Fix the individual montage size and determine the number
montWidth = 25; % in pixel
% find a traction peak
[~,peakFrameForce] = max(curTrack.forceMag(chosenFRange));
peakFrameForce = chosenFRange(peakFrameForce);
if ~isempty(curTrack.forcePeakness) && curTrack.forcePeakness
    peakFrameForce = curTrack.forcePeakFrame;
end
frameFII=round(curTrack.firstIncreseTimeInt/tInterval);
frameFTI=round(curTrack.firstIncreseTimeForce/tInterval);

maxMontageNum = floor(535/montWidth)*2;
numChosenFrames = length(chosenStartFrame:chosenEndFrame);
montInterval = ceil(numChosenFrames/maxMontageNum);
if montInterval>1
    %Check if peakFrameForce is included in the range
    modShift=mod(peakFrameForce-chosenStartFrame,montInterval);
    indiceRange=1+modShift:montInterval:numChosenFrames;
    % Add correct first rise frames
    frameFIIcut = frameFII -chosenStartFrame +1;
    frameFTIcut = frameFTI -chosenStartFrame +1;
    peakFrameCut = peakFrame-chosenStartFrame +1;
    if ~ismember(frameFIIcut,indiceRange)
        indiceRange = [indiceRange frameFIIcut];
        indiceRange = sort(indiceRange);
    end
    if ~ismember(frameFTIcut,indiceRange)
        indiceRange = [indiceRange frameFTIcut];
        indiceRange = sort(indiceRange);
    end
    if ~ismember(peakFrameCut,indiceRange)
        indiceRange = [indiceRange peakFrameCut];
        indiceRange = sort(indiceRange);
    end
elseif montInterval==1
    indiceRange=1:montInterval:numChosenFrames;
end
    
hm1=montage(cropImg,'Size',[2, NaN],'DisplayRange',[minInt maxInt], 'Indices', indiceRange);
line([2 2+scaleBar*500/pixSize],...
    [3 3],'LineWidth',2,'Color',[1,1,1]);
monImgH = floor(hm1.YData(2)/2);
numCols = ceil(length(indiceRange)/2);
monImgW = floor(hm1.XData(2)/numCols);
p=0;
hold on

for ii= indiceRange
    p=p+1;
    iCol=mod(p-1,numCols);
    q=floor((p-1)/numCols);
    txtMont1= text(1+(iCol)*monImgW, (q+1)*monImgH-1,[num2str((ii-1)*tInterval,'% 10.0f') ' s'],'Color',[1,1,1]);
    set(txtMont1,'FontUnits','pixels')
    set(txtMont1,'Fontsize',6)
    set(txtMont1,'horizontalAlignment','right')
    set(txtMont1,'position',[(iCol+1)*monImgW-0.5, (q+1)*monImgH-1.5])
    % Showing tracks
%     sF = chosenStartFrame;
    eF = chosenStartFrame+ii-1;
%     plot(curTrack.xCoord(sF:eF)-bLeft+(iCol)*monImgW+1,curTrack.yCoord(sF:eF)-bBottom+q*monImgH+1,'r', 'LineWidth', 0.5)
    % Show first increase point if it exists
    if eF>=curStartFrame && eF<curStartFrame+montInterval
        markerType = 'yo';
    elseif eF==peakFrame
        markerType = 'wo';
    elseif eF>=curEndFrame && eF<curEndFrame+montInterval
        markerType = 'bo';
    else
        markerType = 'ro';
    end
    if ~isempty(curTrack.forceTransmitting) && curTrack.forceTransmitting && (eF==frameFII)% && eF<frameFII+montInterval)
        markerType = 'go';
    end
    plot(curTrack.xCoord(eF)-bLeft+(iCol)*monImgW+1,curTrack.yCoord(eF)-bBottom+q*monImgH+1,markerType,'MarkerSize',7, 'LineWidth', 0.5)
%         plot(curTrack.xCoord(chosenStartFrame:(chosenStartFrame+ii-1))-bLeft+(iCol)*monImgW+1,curTrack.yCoord(chosenStartFrame:(chosenStartFrame+ii-1))-bBottom+q*monImgH+1,'r', 'LineWidth', 0.5)
%         plot(curTrack.xCoord(chosenStartFrame+ii-1)-bLeft+(iCol)*monImgW+1,curTrack.yCoord(chosenStartFrame+ii-1)-bBottom+q*monImgH+1,'ro','MarkerSize',8, 'LineWidth', 0.5)
end

ax7=subplot('Position',[marginX, 200/figHeight, 540/figWidth-marginX,50/figHeight]);
% tCropMax = max(curTrack.forceMag(chosenStartFrame:chosenEndFrame))*1.1;
% tCropMin = min(curTrack.forceMag(chosenStartFrame:chosenEndFrame));
try
    montage(cropTmap,'DisplayRange',[tCropMin, tCropMax],'Size',[2, NaN], 'Indices', indiceRange)
catch
    montage(cropTmap,'DisplayRange',[tCropMin, 10],'Size',[2, NaN], 'Indices', indiceRange)
end
p=0;
hold on
for ii= indiceRange
    p=p+1;
    iCol=mod(p-1,numCols);
    q=floor((p-1)/numCols);
    txtMont2= text(1+(iCol)*monImgW, (q+1)*monImgH-1,[num2str((ii-1)*tInterval,'% 10.0f') ' s'],'Color','w');
    set(txtMont2,'FontUnits','pixels')
    set(txtMont2,'Fontsize',6)
    set(txtMont2,'horizontalAlignment','right')
    set(txtMont2,'position',[(iCol+1)*monImgW-0.5, (q+1)*monImgH-1.5])
    % Showing tracks
%     sF = chosenStartFrame;
    eF = chosenStartFrame+ii-1;
    if eF>=curStartFrame && eF<curStartFrame+montInterval
        markerType = 'yo';
    elseif eF==peakFrameForce 
        markerType = 'wo';
    elseif eF>=curEndFrame && eF<curEndFrame+montInterval
        markerType = 'bo';
    else
        markerType = 'ro';
    end
    if ~isempty(curTrack.forceTransmitting) && curTrack.forceTransmitting && (eF==frameFTI)% && eF<frameFTI+montInterval)
        markerType = 'go';
    end
%     plot(curTrack.xCoord(sF:eF)-bLeft+(iCol)*monImgW+1,curTrack.yCoord(sF:eF)-bBottom+q*monImgH+1,'w', 'LineWidth', 0.5)
    plot(curTrack.xCoord(eF)-bLeft+(iCol)*monImgW+1,curTrack.yCoord(eF)-bBottom+q*monImgH+1,markerType,'MarkerSize',7, 'LineWidth', 0.5)
end
colormap(ax4,'jet')
colormap(ax5,'jet')
colormap(ax6,'jet')
colormap(ax7,'jet')
% intensity plot
%     subplot('Position',[0,0,0.5,0.33]), plot(chosenStartFrame:chosenEndFrame,curTrack.ampTotal(chosenStartFrame:chosenEndFrame))
%     subplot('Position',[0.5,0,0.5,0.33]), plot(chosenStartFrame:chosenEndFrame,curTrack.forceMag(chosenStartFrame:chosenEndFrame),'r')
% ax 8: intensity time series

ax8=axes('Position',[60/figWidth, 50/figHeight, 200/figWidth-marginX,130/figHeight]);
plot((curStartFrameEE-curStartFrameEE:curEndFrameEE-curStartFrameEE)*tInterval,curTrack.ampTotal(curStartFrameEE:curEndFrameEE),'k'), hold on
plot((curStartFrame-curStartFrameEE:curEndFrame-curStartFrameEE)*tInterval,curTrack.ampTotal(curStartFrame:curEndFrame),'b')
if ~isempty(curTrack.forceTransmitting) && curTrack.forceTransmitting
    plot((frameFII-curStartFrameEE)*tInterval,curTrack.ampTotal(frameFII),'o','MarkerFaceColor','b','MarkerEdgeColor','w')
end
%background level
line([0 (curEndFrameEE-curStartFrameEE)*tInterval],[curTrack.bkgMaxInt curTrack.bkgMaxInt],'linestyle',':','Color','k')
xlabel('Time (s)','FontUnits',genFontUnit,'FontSize',genFontSize); ylabel('Fluorescence intensity (a.u.)','FontUnits',genFontUnit,'FontSize',genFontSize)
set(ax8,'FontUnits',genFontUnit,'FontSize',genFontSize)
% force time series
ax9=axes('Position',[350/figWidth, 50/figHeight, 200/figWidth-marginX,130/figHeight]);
plot((curStartFrameEE-curStartFrameEE:curEndFrameEE-curStartFrameEE)*tInterval,curTrack.forceMag(curStartFrameEE:curEndFrameEE),'k'), hold on
plot((curStartFrame-curStartFrameEE:curEndFrame-curStartFrameEE)*tInterval,curTrack.forceMag(curStartFrame:curEndFrame),'r')
if ~isempty(curTrack.forceTransmitting) && curTrack.forceTransmitting
    plot((frameFTI-curStartFrameEE)*tInterval,curTrack.forceMag(frameFTI),'o','MarkerFaceColor','r','MarkerEdgeColor','w')
end
try
    line([0 (curEndFrameEE-curStartFrameEE)*tInterval],[curTrack.bkgMaxForce curTrack.bkgMaxForce],'LineStyle',':','Color','k')
catch
    curTrack = calculateFirstIncreaseTimeTracks(curTrack,0.5,0.05,tInterval);
    try
        line([0 (curEndFrameEE-curStartFrameEE)*tInterval],[curTrack.bkgMaxForce curTrack.bkgMaxForce],...
            'LineStyle',':','Color','k')
    catch
        disp(' ')
    end
end
xlabel('Time (s)'); ylabel('Traction (Pa)')
set(ax9,'FontUnits',genFontUnit,'FontSize',genFontSize)

% Smoothed line and maximum
splineParam = 0.01;
d = curTrack.ampTotal(curStartFrameEE:curEndFrameEE);
tRange = curTrack.iFrame(curStartFrameEE:curEndFrameEE);
%     numNan = find(isnan(d),1,'last');
%     tRange(isnan(d)) = [];
%     d(isnan(d)) = [];
axes(ax8)
sd_spline= csaps(tRange,d,splineParam);
sd=ppval(sd_spline,tRange);
plot((curStartFrameEE-curStartFrameEE:curEndFrameEE-curStartFrameEE)*tInterval,sd,'Color',[0/255 0/255 0/255],'Linewidth',2)
plot((curStartFrame-curStartFrameEE:curEndFrame-curStartFrameEE)*tInterval,...
    sd(curStartFrame-curStartFrameEE+1:curEndFrame-curStartFrameEE+1),'Color',[84/255 84/255 255/255],'Linewidth',2)
% Plot maximum point
if ~isempty(curTrack.intenPeakness) && curTrack.intenPeakness
    plot((curTrack.intenPeakFrame-curStartFrameEE)*tInterval,sd(curTrack.intenPeakFrame-curStartFrameEE+1),'o',...
        'MarkerFaceColor','w','MarkerEdgeColor',[84/255 84/255 255/255])
end


axes(ax9)
curForce = curTrack.forceMag(curStartFrameEE:curEndFrameEE);
sCurForce_spline= csaps(tRange,curForce,splineParam);
sCurForce=ppval(sCurForce_spline,tRange);
plot((curStartFrameEE-curStartFrameEE:curEndFrameEE-curStartFrameEE)*tInterval,sCurForce,'Color',[0/255 0/255 0/255],'Linewidth',2)
plot((curStartFrame-curStartFrameEE:curEndFrame-curStartFrameEE)*tInterval,...
    sCurForce(curStartFrame-curStartFrameEE+1:curEndFrame-curStartFrameEE+1),'Color',[229/255 84/255 84/255],'Linewidth',2)
% Plot maximum point
if ~isempty(curTrack.forcePeakness) && curTrack.forcePeakness
    plot((curTrack.forcePeakFrame-curStartFrameEE)*tInterval,sCurForce(curTrack.forcePeakFrame-curStartFrameEE+1),'o',...
        'MarkerFaceColor','w','MarkerEdgeColor',[229/255 84/255 84/255])
end

% force map colorbar
ax10=axes('Position',[3*marginX+3*175/figWidth, 260/figHeight+marginY, 60/figWidth,145/figHeight]);
caxis([tmin tmax]); axis off
%     if isempty(hc)
axis tight
hcb1 = colorbar('West');
hcb1.Position = [ax10.Position(1)+10/figWidth ax10.Position(2)*1.01 ax10.Position(3)*0.15 ax10.Position(4)*0.9];
set(ax10,'FontUnits',genFontUnit,'FontSize',genFontSize)
set(hcb1,'YAxisLocation','right')
set(get(hcb1,'xlabel'),'String','Traction (Pa)')

% force colorbar mini for montage
ax11=axes('Position',[540/figWidth, 203/figHeight, 60/figWidth, 50/figHeight]);
axis tight
caxis([tCropMin, tCropMax]), axis off
%     if isempty(hc)
hcb2 = colorbar('West');
hcb2.Position = [ax11.Position(1)+10/figWidth ax11.Position(2)*1.01 ax11.Position(3)*0.15 ax11.Position(4)*0.9];
set(ax11,'FontUnits',genFontUnit,'FontSize',genFontSize)
set(hcb2,'YAxisLocation','right')
set(get(hcb2,'xlabel'),'String','Traction (Pa)')
colormap(ax10,'jet')
colormap(ax11,'jet')

if exist('gPath','var')
    print(h2,strcat(gPath,'/track',num2str(IDtoInspect),additionalName,'.eps'),'-depsc2')
    savefig(h2,strcat(gPath,'/track',num2str(IDtoInspect),additionalName,'.fig'))
end
if exist('IDtoInspect','var')
    disp(['This ' num2str(IDtoInspect) ' th track has ' num2str(curTrack.firstIncreseTimeIntAgainstForce,'% 10.1f') ' sec of initial time lag of force after fluorescence signal.'])
end
disp('Marker: yellow (track start), green (initial rise), blue(track end), white (peak)')
end
function [h2] = showSingleAdhesionTrackSummaryRateConstFitting(MD,curTrack,imgMap,tMap,imgMap2,IDtoInspect, gPath,additionalName)
% h2 = showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,IDtoInspect, gPath,additionalName)
% This function shows big picture, montage, and time series of fluorescence
% intensity and traction.
% Sangyoon Han, Jun 2015
imgWidth = size(imgMap,2);
imgHeight = size(imgMap,1);

pixSize = MD.pixelSize_; % nm/pixel
tInterval = MD.timeInterval_; % time interval in sec
scaleBar = 1; %micron

curEndFrame=curTrack.endingFrameExtra;
curStartFrame = curTrack.startingFrameExtra;
curEndFrameEE=curTrack.endingFrameExtraExtra;
curStartFrameEE = max(1,curStartFrame-5); %This needs to be updated in calculatePeakTimeLagFromTracks %curTrack.startingFrameExtraExtra;
curFrameRange= curStartFrameEE:curEndFrameEE;
chosenStartFrame = curStartFrameEE;
chosenEndFrame = curEndFrameEE;

chosenFRange = curFrameRange;

splineParam = 0.01;
try
    d = curTrack.amp(curStartFrameEE:curEndFrameEE);
catch
%     curTrack=readIntensityFromTracks(curTrack,imgMap,1,'extraLength',120,'movieData',MD,'reTrack',false);
%     curTrack=readIntensityFromTracks(curTrack,tMap,2,'extraLength',120,'movieData',MD);
%     if ~isempty(imgMap2)
%         curTrack=readIntensityFromTracks(curTrack,imgMap2,5,'extraLength',120,'movieData',MD);
%     end   
    d = curTrack.ampTotal(curStartFrame:curEndFrame);
end    
tRange = curTrack.iFrame(curStartFrameEE:curEndFrameEE);
sd_spline= csaps(tRange,d,splineParam);
sd=ppval(sd_spline,tRange);

[~,peakFrame] = max(sd); %max(curTrack.ampTotal(chosenFRange));
peakFrame = chosenFRange(peakFrame);

if isfield(curTrack,'intenPeakness') && ~isempty(curTrack.intenPeakness) && curTrack.intenPeakness
    peakFrame = curTrack.intenPeakFrame;
end
if peakFrame >= ((chosenFRange(end)-chosenFRange(1))*0.8 + chosenFRange(1))
    middleFrame = round((chosenFRange(1)+chosenFRange(end))/2);
else
    middleFrame = peakFrame;
end

maxX = nanmax(curTrack.xCoord)-nanmin(curTrack.xCoord);
maxY =  nanmax(curTrack.yCoord)-nanmin(curTrack.yCoord);
r_pix = ceil(max(max(maxX,maxY)/2+5,10));
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
h2=figure('Name', strcat('trackID',num2str(IDtoInspect)));
% Depending on input, figure might have to be bigger
figWidth = 600+200; % addtinal 200 is for new graph place.
figHeight = 675;
if ~isempty(imgMap2)
    figHeight = figHeight + 175 + 5 + 50 + 5;
end
%     set(h2,'Position',[100,10,figWidth,figHeight]),title(['ID:' num2str(IDtoInspect)])
set(h2,'Position',[30,100,figWidth,figHeight])
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
midStartFrame = chosenStartFrame; %ceil((chosenStartFrame+peakFrame)/2);
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
% ax2: peak point or middle point
ax2=subplot('Position',[2*marginX+175/figWidth, 490/figHeight, 175/figWidth, 175/figHeight]);
imshow(imgMap(:,:,middleFrame),[minInt maxInt]), hold on
rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
set(ax2,'XLim',imgBigXLim,'YLim',imgBigYLim)
line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
    [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
peakFrameAdjusted = middleFrame - chosenStartFrame; % in frames
txtAx2= text(bRightBig-40, bBottomBig+12,[num2str(peakFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
set(txtAx2,'horizontalAlignment','right')
set(txtAx2,'position',[bRightBig-mgScale,bBottomBig+12])
set(findobj(ax2,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)
% ax3: end point
% it would be better to show frame that actually shows adhesion, so we
% will use frame that is mid-point between initial and peak
midEndFrame = chosenEndFrame; %ceil((chosenEndFrame+peakFrame)/2);
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
montWidth = 25; % in pixel
if ~isempty(tMap)
    tMap2 = tMap(:,:,[midStartFrame,middleFrame,midEndFrame]);
    tMapROI = tMap2(bBottomBig:bTopBig,bLeftBig:bRightBig,:);
    % tmax= min(tCropMax*3,max(tMapROI(:))*0.9);
    % tmin= tCropMin*0.1;
    tmax= max(tMapROI(:))*0.9;
    tmin= min(tMapROI(:));

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
    cropTmap = reshape(cropTmap,size(cropTmap,1),size(cropTmap,2),1,size(cropTmap,3));
    % Fix the individual montage size and determine the number
    % tCropMax = max(curTrack.forceMag(chosenStartFrame:chosenEndFrame))*1.1;
    % tCropMin = min(curTrack.forceMag(chosenStartFrame:chosenEndFrame));
    tCropMax = max(cropTmap(:))*0.9;
    tCropMin = min(cropTmap(:));
    % find a traction peak
    [~,peakFrameForce] = max(curTrack.forceMag(chosenFRange));
    peakFrameForce = chosenFRange(peakFrameForce);
    if isfield(curTrack,'forcePeakness') && curTrack.forcePeakness
        peakFrameForce = curTrack.forcePeakFrame;
    end
end
% make 4D array to use montage
cropImg = reshape(cropImg,size(cropImg,1),size(cropImg,2),1,size(cropImg,3));
axes('Position',[marginX, 430/figHeight+marginY, 540/figWidth-marginX,50/figHeight]);
    
maxMontageNum = floor(535/montWidth/2); %*2;
numChosenFrames = length(chosenStartFrame:chosenEndFrame);
montInterval = ceil(numChosenFrames/maxMontageNum);
if montInterval>1 && ~isempty(tMap)
    %Check if peakFrameForce is included in the range
    modShift=mod(peakFrameForce-chosenStartFrame,montInterval);
    indiceRange=1+modShift:montInterval:numChosenFrames;
else
    indiceRange=1:montInterval:numChosenFrames;
end
    
hm1=montage(cropImg,'Size',[1, NaN],'DisplayRange',[minInt maxInt], 'Indices', uint16(indiceRange),'ThumbnailSize',[2*r_pix+1 2*r_pix+1]);
line([2 2+scaleBar*500/pixSize],...
    [3 3],'LineWidth',2,'Color',[1,1,1]);
monImgH = floor(hm1.YData(2)); %/2);
numCols = ceil(length(indiceRange)); %/2);
monImgW = floor(hm1.XData(2)/numCols);
p=0;
hold on

for ii= indiceRange
    p=p+1;
    iCol=mod(p-1,numCols);
    q=floor((p-1)/numCols);
    txtMont1= text(1+(iCol)*monImgW, (q+1)*monImgH-4,[num2str((ii-1)*tInterval,'% 10.0f') ' s'],'Color',[1,1,1]);
    set(txtMont1,'FontUnits','pixels')
    set(txtMont1,'Fontsize',7)
    set(txtMont1,'horizontalAlignment','right')
    set(txtMont1,'position',[(iCol+1)*monImgW-0.5, (q+1)*monImgH-1.5])
    % Showing tracks
    eF = chosenStartFrame+ii-1;
    % Show first increase point if it exists
    if eF>=curStartFrame && eF<curStartFrame+montInterval
        markerType = 'yo'; bdType='y';
    elseif eF==peakFrame
        markerType = 'wo'; bdType='w';
    elseif eF>=curEndFrame && eF<curEndFrame+montInterval
        markerType = 'bo'; bdType='b';
    else
        markerType = 'ro'; bdType='r';
    end
    if isfield(curTrack,'forceTransmitting') && curTrack1.forceTransmitting && (eF==frameFII)% && eF<frameFII+montInterval)
        markerType = 'go'; bdType='g';
    end
    plot(curTrack.xCoord(eF)-bLeft+(iCol)*monImgW+1,curTrack.yCoord(eF)-bBottom+q*monImgH+1,markerType,'MarkerSize',14, 'LineWidth', 0.5)
end

% Montage plot for tMap
if ~isempty(tMap)
    ax7=subplot('Position',[marginX, 200/figHeight, 540/figWidth-marginX,50/figHeight]);
    % tCropMax = max(curTrack.forceMag(chosenStartFrame:chosenEndFrame))*1.1;
    % tCropMin = min(curTrack.forceMag(chosenStartFrame:chosenEndFrame));
    try
        montage(cropTmap,'DisplayRange',[tCropMin, tCropMax],'Size',[1, NaN], 'Indices', indiceRange,'ThumbnailSize',[2*r_pix+1 2*r_pix+1]);
    catch
        montage(cropTmap,'DisplayRange',[tCropMin, 10],'Size',[1, NaN], 'Indices', indiceRange,'ThumbnailSize',[2*r_pix+1 2*r_pix+1]);
    end
    p=0;
    hold on
    for ii= indiceRange
        p=p+1;
        iCol=mod(p-1,numCols);
        q=floor((p-1)/numCols);
        txtMont2= text(1+(iCol)*monImgW, (q+1)*monImgH-4,[num2str((ii-1)*tInterval,'% 10.0f') ' s'],'Color','w');
        set(txtMont2,'FontUnits','pixels')
        set(txtMont2,'Fontsize',7)
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
        if isfield(curTrack,'forceTransmitting') && curTrack1.forceTransmitting && (eF==frameFTI)% && eF<frameFTI+montInterval)
            markerType = 'go';
        end
        plot(curTrack.xCoord(eF)-bLeft+(iCol)*monImgW+1,curTrack.yCoord(eF)-bBottom+q*monImgH+1,markerType,'MarkerSize',14, 'LineWidth', 0.5)
    end
    colormap(ax4,'jet')
    colormap(ax5,'jet')
    colormap(ax6,'jet')
    colormap(ax7,'jet')
end
% intensity plot
%     subplot('Position',[0,0,0.5,0.33]), plot(chosenStartFrame:chosenEndFrame,curTrack.ampTotal(chosenStartFrame:chosenEndFrame))
%     subplot('Position',[0.5,0,0.5,0.33]), plot(chosenStartFrame:chosenEndFrame,curTrack.forceMag(chosenStartFrame:chosenEndFrame),'r')
% ax 8: intensity time series

if ~isempty(imgMap2)
    ax8=axes('Position',[4*marginX+(3*175+60+30)/figWidth, (430+80+40)/figHeight, 155/figWidth,80/figHeight]);
else
    ax8=axes('Position',[50/figWidth, 50/figHeight, 150/figWidth-marginX,130/figHeight]);
end
plot((curStartFrame-curStartFrameEE:curEndFrame-curStartFrameEE)*tInterval,curTrack.amp(curStartFrame:curEndFrame),'k'), hold on    
xlabel('Time (s)','FontUnits',genFontUnit,'FontSize',genFontSize); ylabel('Fluorescence intensity (a.u.)','FontUnits',genFontUnit,'FontSize',genFontSize)
set(ax8,'FontUnits',genFontUnit,'FontSize',genFontSize)

if ~isempty(imgMap2)
    axA=axes('Position',[4*marginX+(3*175+60+30)/figWidth, (330+80+40)/figHeight, 155/figWidth,80/figHeight]);
else
    axA=axes('Position',[50/figWidth, 200/figHeight, 150/figWidth-marginX,130/figHeight]);
end
tIntervalMin = MD.timeInterval_/60;
[assemRate,bestModelAssem,bestSummaryAssem,yEntireAssem] = getAssemRateFromTracks(curTrack,tIntervalMin);
% bestModelAssem.plot
if ~isempty(bestModelAssem)
    plot(bestModelAssem.Variables.x1*60,bestModelAssem.Variables.y, 'bx-.')
    hold on
    plot(bestModelAssem.Variables.x1*60,bestModelAssem.Fitted, 'r-')

    xlabel('Time (s)','FontUnits',genFontUnit,'FontSize',genFontSize); 
    ylabel('Log[I/Io]','FontUnits',genFontUnit,'FontSize',genFontSize)
    text(bestModelAssem.Variables.x1(1)*60, (bestModelAssem.Fitted(1)+2*bestModelAssem.Fitted(end))/3, {['R^2_{adj}=' num2str(bestSummaryAssem.adjRsquared) ];['K_{assem}=' num2str(assemRate)]})
else
    text(0,0.5,{'Failed assembly rate'; 'estimation'})
end
set(axA,'FontUnits',genFontUnit,'FontSize',genFontSize)


if ~isempty(imgMap2)
    axB=axes('Position',[4*marginX+(3*175+60+30)/figWidth, (240+80+40)/figHeight, 155/figWidth,70/figHeight]);
else
    axB=axes('Position',[50/figWidth, 100/figHeight, 150/figWidth-marginX,80/figHeight]);
end
tIntervalMin = MD.timeInterval_/60;
[disassemRate,bestModelDis,bestSummaryDis,yEntireDis] = getDisassemRateFromTracks(curTrack,tIntervalMin);
% bestModelAssem.plot
if ~isempty(bestModelDis)
    plot(bestModelDis.Variables.x1*60,bestModelDis.Variables.y, 'bx-.')
    hold on
    plot(bestModelDis.Variables.x1*60,bestModelDis.Fitted, 'r-')
    
    xlabel('Time (s)','FontUnits',genFontUnit,'FontSize',genFontSize); 
    ylabel('Log[Io/I]','FontUnits',genFontUnit,'FontSize',genFontSize)
    text(bestModelDis.Variables.x1(2)*60, (bestModelDis.Fitted(1)+2*bestModelDis.Fitted(end))/3, {['R^2_{adj}=' num2str(bestSummaryDis.adjRsquared) ];['K_{disassem}=' num2str(disassemRate)]})
else
    text(0,0.5,{'Failed disassembly rate'; 'estimation'})
end
set(axB,'FontUnits',genFontUnit,'FontSize',genFontSize)

% force time series
if ~isempty(tMap)
    if ~isempty(imgMap2)
        ax9=axes('Position',[4*marginX+(3*175+60+30)/figWidth, (120+80+70)/figHeight, 155/figWidth,70/figHeight]);
    else
        ax9=axes('Position',[250/figWidth, 50/figHeight, 150/figWidth-marginX,90/figHeight]);
    end
    plot((curStartFrameEE-curStartFrameEE:curEndFrameEE-curStartFrameEE)*tInterval,curTrack.forceMag(curStartFrameEE:curEndFrameEE),'k'), hold on
    plot((curStartFrame-curStartFrameEE:curEndFrame-curStartFrameEE)*tInterval,curTrack.forceMag(curStartFrame:curEndFrame),'r')
    xlabel('Time (s)'); ylabel('Traction (Pa)')
    set(ax9,'FontUnits',genFontUnit,'FontSize',genFontSize)

    % Cross variance plot
    if ~isempty(imgMap2)
        ax12=axes('Position',[4*marginX+(3*175+60+30)/figWidth, (200+30)/figHeight, 155/figWidth,80/figHeight]);
    else
        ax12=axes('Position',[450/figWidth, 50/figHeight, 150/figWidth-marginX,130/figHeight]);
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
end
%% Third channel plotting
if ~isempty(imgMap2)
    % ax13: initial point for the third channel
    ax13=axes('Position',[marginX, (490+175+5+50+5)/figHeight, 175/figWidth, 175/figHeight]);
    cropImg2 = imgMap2(bBottom:bTop,bLeft:bRight,chosenStartFrame:chosenEndFrame);
    cropImgSmaller = imgMap2(bBottom+7:bTop-7,bLeft+7:bRight-7,chosenStartFrame:chosenEndFrame);
    maxInt2 = max(cropImgSmaller(:));
    minInt2 = min(cropImgSmaller(:));
    imshow(imgMap2(:,:,midStartFrame),[minInt2 maxInt2]), hold on
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax13,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    txtAx1= text(bRightBig-40, bBottomBig+12,[num2str(midStartFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
    set(txtAx1,'horizontalAlignment','right')
    set(txtAx1,'position',[bRightBig-mgScale,bBottomBig+12])
    set(findobj(ax13,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)

    % ax14: peak point or middle point
    ax14=subplot('Position',[2*marginX+175/figWidth, (490+175+5+50+5)/figHeight, 175/figWidth, 175/figHeight]);
    imshow(imgMap2(:,:,middleFrame),[minInt2 maxInt2]), hold on
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax14,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    peakFrameAdjusted = middleFrame - chosenStartFrame; % in frames
    txtAx2= text(bRightBig-40, bBottomBig+12,[num2str(peakFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
    set(txtAx2,'horizontalAlignment','right')
    set(txtAx2,'position',[bRightBig-mgScale,bBottomBig+12])
    set(findobj(ax14,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)

    % ax15: end point
    ax15=subplot('Position',[3*marginX+2*175/figWidth, (490+175+5+50+5)/figHeight, 175/figWidth, 175/figHeight]);
    imshow(imgMap2(:,:,midEndFrame),[minInt2 maxInt2]), hold on
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax15,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    midEndFrameAdjusted = midEndFrame - chosenStartFrame; % in frames
    txtAx3= text(bRightBig-40, bBottomBig+12,[num2str(midEndFrameAdjusted*tInterval,'% 10.0f') ' sec'],'Color',[1,1,1]);
    set(txtAx3,'horizontalAlignment','right')
    set(txtAx3,'position',[bRightBig-mgScale,bBottomBig+12])
    set(findobj(ax15,'Type','text'),'FontUnits',genFontUnit,'FontSize',genFontSize)

    %Montage
    axes('Position',[marginX, (430+50+10+175)/figHeight+marginY, 540/figWidth-marginX,50/figHeight]);
    cropImg2 = reshape(cropImg2,size(cropImg2,1),size(cropImg2,2),1,size(cropImg2,3));
    hm3=montage(cropImg2,'Size',[1, NaN],'DisplayRange',[minInt2 maxInt2], 'Indices', uint16(indiceRange),'ThumbnailSize',[2*r_pix+1 2*r_pix+1]);
    hold on
    p=0;
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
        eF = chosenStartFrame+ii-1;
        % Show first increase point if it exists
        if eF>=curStartFrame && eF<curStartFrame+montInterval
            markerType = 'yo'; bdType='y';
        elseif eF==peakFrame
            markerType = 'wo'; bdType='w';
        elseif eF>=curEndFrame && eF<curEndFrame+montInterval
            markerType = 'bo'; bdType='b';
        else
            markerType = 'ro'; bdType='r';
        end
        plot(curTrack.xCoord(eF)-bLeft+(iCol)*monImgW+1,curTrack.yCoord(eF)-bBottom+q*monImgH+1,markerType,'MarkerSize',14, 'LineWidth', 0.5)
    end

    % Third channel time series

    ax16=axes('Position',[4*marginX+(3*175+60+30)/figWidth, (490+175+80+70)/figHeight, 155/figWidth,80/figHeight]);
%     plot((curStartFrameEE-curStartFrameEE:curEndFrameEE-curStartFrameEE)*tInterval,curTrack.amp2(curStartFrameEE:curEndFrameEE),'k'), hold on

    try
        plot((curStartFrame-curStartFrameEE:curEndFrame-curStartFrameEE)*tInterval,curTrack.amp2(curStartFrame:curEndFrame),'m')
    catch
        plot((curStartFrame-curStartFrameEE:curEndFrame-curStartFrameEE)*tInterval,curTrack.ampTotal2(curStartFrame:curEndFrame),'m')
    end
    
    xlabel('Time (s)'); ylabel('amp 2 (A.U.)')
    set(ax16,'FontUnits',genFontUnit,'FontSize',genFontSize)

    % Assembly for amp2
    ax17=axes('Position',[4*marginX+(3*175+60+30)/figWidth, (410+175+80+70)/figHeight, 155/figWidth,60/figHeight]);
    [assemRate2,bestModelAssem2,bestSummaryAssem2] = getAssemRateFromTracks(curTrack,tIntervalMin,'amp2');
    % bestModelAssem.plot
    if ~isempty(bestModelAssem2)
        plot(bestModelAssem2.Variables.x1*60,bestModelAssem2.Variables.y, 'bx-.')
        hold on
        plot(bestModelAssem2.Variables.x1*60,bestModelAssem2.Fitted, 'r-')

        xlabel('Time (s)','FontUnits',genFontUnit,'FontSize',genFontSize); 
        ylabel('Log[I/Io]','FontUnits',genFontUnit,'FontSize',genFontSize)
        text(bestModelAssem2.Variables.x1(1)*60, (bestModelAssem2.Fitted(1)+2*bestModelAssem2.Fitted(end))/3, {['R^2_{adj}=' num2str(bestSummaryAssem2.adjRsquared) ];['K_{assem}=' num2str(assemRate2)]})
    else
        text(0,0.5,{'Failed assembly rate'; 'estimation'})
    end
    set(ax17,'FontUnits',genFontUnit,'FontSize',genFontSize)

    % Disassembly rate
    ax18=axes('Position',[4*marginX+(3*175+60+30)/figWidth, (320+175+80+70)/figHeight, 155/figWidth,60/figHeight]);
    [disassemRate2,bestModelDis2,bestSummaryDis2,~] = getDisassemRateFromTracks(curTrack,tIntervalMin,'amp2');
    % bestModelAssem.plot
    if ~isempty(bestModelDis2)
        plot(bestModelDis2.Variables.x1*60,bestModelDis2.Variables.y, 'bx-.')
        hold on
        plot(bestModelDis2.Variables.x1*60,bestModelDis2.Fitted, 'r-')

        xlabel('Time (s)','FontUnits',genFontUnit,'FontSize',genFontSize); 
        ylabel('Log[Io/I]','FontUnits',genFontUnit,'FontSize',genFontSize)
        text(bestModelDis2.Variables.x1(2)*60, (bestModelDis2.Fitted(1)+2*bestModelDis2.Fitted(end))/3, {['R^2_{adj}=' num2str(bestSummaryDis2.adjRsquared) ];['K_{disassem}=' num2str(disassemRate2)]})
    else
        text(0,0.5,{'Failed disassembly rate'; 'estimation'})
    end
    set(ax18,'FontUnits',genFontUnit,'FontSize',genFontSize)
    
else
    timeLagMasterAgainstMainSlave = NaN;
end
%% saving

if exist('gPath','var')
    print(h2,strcat(gPath,'/trackID',num2str(IDtoInspect),'rateConstant.eps'),'-depsc2')
    print(h2,strcat(gPath,'/trackID',num2str(IDtoInspect),'rateConstant.png'),'-dpng')
    savefig(h2,strcat(gPath,'/trackID',num2str(IDtoInspect),'rateConstant.fig'))
end
disp('Marker: yellow (track start), green (initial rise), blue(track end), white (peak)')
end
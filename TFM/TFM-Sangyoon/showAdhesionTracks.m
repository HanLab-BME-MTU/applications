function [IDs]=showAdhesionTracks(pathForColocalization,idList,varargin)
% this function reads from Colocalization folder and shows the specified
% tracks on selected Channel interactively.
%% input reading
ip =inputParser;
ip.addRequired('pathForColocalization',@ischar)
ip.addOptional('idList','all',@(x)islogical(x)||ischar(x)||isempty(x))
ip.addParamValue('numChan',2,@isscalar); % selcted track ids
ip.addParamValue('tracksNA',[],@isstruct); % selcted track ids
ip.parse(pathForColocalization,idList,varargin{:});
pathForColocalization=ip.Results.pathForColocalization;
idList=ip.Results.idList;
tracksNA=ip.Results.tracksNA;
numChan=ip.Results.numChan;

%% Load processed data
disp('Loading raw files ...')
tic
if isempty(tracksNA)
    tracksNA = load([pathForColocalization filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNA = tracksNA.tracksNA;
end
if numChan==1
    imgMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
    imgMap = imgMap.tMap;
elseif numChan==2
    imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
    imgMap = imgMap.paxImgStack;
    tMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
    tMap = tMap.tMap;
end
outputPath = [pathForColocalization filesep 'trackAnalysis'];
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end
numFrames = size(imgMap,3);
if ischar(idList) && strcmp(idList,'all')
    startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA)));
    endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA)));
    idxIDList = 1:numel(tracksNA);
else
    tracksNA=tracksNA(idList);
    idxIDList = find(idList);
    startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA)));
    endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA)));
%     startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA(idList))));
%     endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA(idList))));
end
% movieData to find out pixel size
coloPath = fileparts(pathForColocalization);
MDPath = fileparts(coloPath);
MDfilePath = [MDPath filesep 'movieData.mat'];
MD = load(MDfilePath,'MD');
MD = MD.MD;
pixSize = MD.pixelSize_; % nm/pixel
tInterval = MD.timeInterval_; % time interval in sec
scaleBar = 1; %micron
toc
%% show interactive movie interface
curNumFrames = endFrame-startFrame+1;
hFig = figure('Position',[100 100 720 750],'Units','normalized','DeleteFcn',@windlowClose);

handles.axes1 = axes('Units','normalized','Position',[0 0.05 1 0.95]);

%// Create slider and listener object for smooth visualization
handles.SliderFrame = uicontrol('Style','slider','Units','normalized','Position',[0 0 0.6 0.05],...
    'Min',startFrame,'Max',endFrame,'Value',startFrame,'SliderStep',[1/curNumFrames 2/curNumFrames],'Callback',@XSliderCallback);
handles.SliderxListener = handles.SliderFrame.addlistener('Value','PreSet',@(s,e) XListenerCallBack);
% handles.SliderxListener = addlistener(handles.SliderFrame,'Value','ContinuousValueChange',@(s,e) XListenerCallBack);

% handles.Text1 = uicontrol('Style','Text','Position',[180 420 60 30],'String','Current frame');
handles.Edit1 = uicontrol('Style','Edit','Units','normalized','Position',[0.6 0 0.05 0.05],'String',num2str(startFrame));
handles.Text2 = uicontrol('Style','Text','Units','normalized','Position',[0.65 0 0.15 0.05],'String',{'Adhesion ID to inspect'});
handles.Edit2 = uicontrol('Style','Edit','Units','normalized','Position',[0.8 0 0.06 0.05],'String',num2str(startFrame));
handles.PushB2 = uicontrol('Style','pushbutton','Units','normalized','Position',[0.86 0 0.14 0.05],'String','Inspect','Callback',@pushInspectAdhesion);

%// Use setappdata to store the image stack and in callbacks, use getappdata to retrieve it and use it. Check the docs for the calling syntax.
setappdata(hFig,'MyMatrix',imgMap); %// You could use %//setappdata(0,'MyMatrix',MyMatrix) to store in the base workspace. 
setappdata(hFig,'tracksNA',tracksNA); 
%// Display 1st frame
imshow(imgMap(:,:,startFrame),[]), hold on
plot(arrayfun(@(x) x.xCoord(startFrame),tracksNA),arrayfun(@(x) x.yCoord(startFrame),tracksNA),'ro')
idAdhLogic = arrayfun(@(x) ~isempty(x.adhBoundary),tracksNA);
idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{startFrame}),tracksNA(idAdhLogic));
idAdh = find(idAdhLogic);
idAdhCur = idAdh(idAdhCur);
arrayfun(@(x) plot(x.adhBoundary{startFrame}(:,1),x.adhBoundary{startFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idAdhCur))
xmat = cell2mat(arrayfun(@(x) x.xCoord(1:startFrame),tracksNA,'UniformOutput',false));
ymat = cell2mat(arrayfun(@(x) x.yCoord(1:startFrame),tracksNA,'UniformOutput',false));
if size(xmat,2)==1
    plot(xmat',ymat','r.')
else
    plot(xmat',ymat','r')
end
hold off
% Supporting data cursor mode to identify an ID of NA track of interest.
dcm_obj = datacursormode(hFig);
set(dcm_obj,'UpdateFcn',@myupdateDC)
imgWidth = size(imgMap,2);
imgHeight = size(imgMap,1);
selectedID = [];
IDs=[];
%// IMPORTANT. Update handles structure.
guidata(hFig,handles);
waitfor(hFig)
if isempty(IDs)
    disp('No track was selected by data cursor. No ID is returend...')
else
    disp(['Selected track is ' num2str(IDs) '.'])
end
function pushInspectAdhesion(~,~)
    IDtoInspect=get(handles.Edit2,'String');
    IDtoInspect = str2double(IDtoInspect);
    % show intensity profile and ask user to pick the starting frame and
    % ending frame
    % read first: 
    IDs=[IDs idxIDList(IDtoInspect)];
    curTrack = readIntensityFromTracks(tracksNA(IDtoInspect),imgMap,1,'extraLength',350);
    curTrack = readIntensityFromTracks(curTrack,tMap,2,'extraLength',350);
    % then use startingFrameExtra and endingFrameExtra to plot intensity
    % time series
    h=figure;
    curEndFrame=curTrack.endingFrameExtra;
    curStartFrame = curTrack.startingFrameExtra;
    curFrameRange= curStartFrame:curEndFrame;
    plot(curFrameRange,curTrack.ampTotal(curFrameRange),'b'), hold on
    plot(curTrack.startingFrame:curTrack.endingFrame,curTrack.ampTotal(curTrack.startingFrame:curTrack.endingFrame),'k')
    set(h,'Position',[500,300,400,200]),title(['ID:' num2str(IDtoInspect)])
    disp('Use data cursor to find out which frame you want to set as starting frame and ending frame')
    chosenStartFrame = input('Chosen starting frame number : ');
    chosenEndFrame = input('Chosen ending frame number : ');
    chosenFRange = chosenStartFrame:chosenEndFrame;
    [~,peakFrame] = max(curTrack.ampTotal(chosenFRange));
    peakFrame = chosenFRange(peakFrame);
    close(h)
    % ROI
%     r_pix = 5; % half the width of the ROI
    % this r_pix was shown to be too small if the maturing adhesion slides.
    % r_pix should include all of the track trace ...
    maxX = nanmax(curTrack.xCoord)-nanmin(curTrack.xCoord);
    maxY =  nanmax(curTrack.yCoord)-nanmin(curTrack.yCoord);
    r_pix = ceil(max(max(maxX,maxY)/2+2,5));
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
    set(h2,'Position',[100,10,figWidth,figHeight])
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
    imshow(imgMap(:,:,midStartFrame),[]), hold on
    % have to show a box indicating a region of interest
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax1,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    midStartFrameAdjusted = midStartFrame - chosenStartFrame; % in frames
    txtAx1= text(bRightBig-40, bBottomBig+12,[num2str(midStartFrameAdjusted*tInterval) ' sec'],'Color',[1,1,1]);
    set(txtAx1,'horizontalAlignment','right')
    set(txtAx1,'position',[bRightBig-mgScale,bBottomBig+12])
    % ax2: peak point
    ax2=subplot('Position',[2*marginX+175/figWidth, 490/figHeight, 175/figWidth, 175/figHeight]);
    imshow(imgMap(:,:,peakFrame),[]), hold on
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax2,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    peakFrameAdjusted = peakFrame - chosenStartFrame; % in frames
    txtAx2= text(bRightBig-40, bBottomBig+12,[num2str(peakFrameAdjusted*tInterval) ' sec'],'Color',[1,1,1]);
    set(txtAx2,'horizontalAlignment','right')
    set(txtAx2,'position',[bRightBig-mgScale,bBottomBig+12])
    % ax3: end point
    % it would be better to show frame that actually shows adhesion, so we
    % will use frame that is mid-point between initial and peak
    midEndFrame = ceil((chosenEndFrame+peakFrame)/2);
    ax3=subplot('Position',[3*marginX+2*175/figWidth, 490/figHeight, 175/figWidth, 175/figHeight]);
    imshow(imgMap(:,:,midEndFrame),[]), hold on
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax3,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    midEndFrameAdjusted = midEndFrame - chosenStartFrame; % in frames
    txtAx3= text(bRightBig-40, bBottomBig+12,[num2str(midEndFrameAdjusted*tInterval) ' sec'],'Color',[1,1,1]);
    set(txtAx3,'horizontalAlignment','right')
    set(txtAx3,'position',[bRightBig-mgScale,bBottomBig+12])
    % ax4: initial point with traction
%     tMap2 = tMap(:,:,[chosenStartFrame,peakFrame,chosenEndFrame]);
    tMap2 = tMap(:,:,[midStartFrame,peakFrame,midEndFrame]);
    tMapROI = tMap2(bBottomBig:bTopBig,bLeftBig:bRightBig,:);
    tmax= max(tMapROI(:))*0.9;
    tmin= min(tMapROI(:));
    
    ax4=subplot('Position',[marginX, 250/figHeight+marginY, 175/figWidth,175/figHeight]);
    imshow(tMap2(:,:,1),[tmin tmax]), hold on
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax4,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    txtAx4= text(bRightBig-40, bBottomBig+12,[num2str(midStartFrameAdjusted*tInterval) ' sec'],'Color',[1,1,1]);
    set(txtAx4,'horizontalAlignment','right')
    set(txtAx4,'position',[bRightBig-mgScale,bBottomBig+12])
    % ax5: fluorescence peak point with traction
    ax5=subplot('Position',[2*marginX+175/figWidth, 250/figHeight+marginY, 175/figWidth,175/figHeight]);
    imshow(tMap2(:,:,2),[tmin tmax]), hold on
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax5,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    txtAx5= text(bRightBig-40, bBottomBig+12,[num2str(peakFrameAdjusted*tInterval) ' sec'],'Color',[1,1,1]);
    set(txtAx5,'horizontalAlignment','right')
    set(txtAx5,'position',[bRightBig-mgScale,bBottomBig+12])
    % ax6: end point with traction
    ax6=subplot('Position',[3*marginX+2*175/figWidth, 250/figHeight+marginY, 175/figWidth,175/figHeight]);
    imshow(tMap2(:,:,3),[tmin tmax]), hold on
    rectangle('Position',[bLeft,bBottom,(bRight-bLeft+1),(bTop-bBottom+1)],'EdgeColor','y'); hold off
    set(ax6,'XLim',imgBigXLim,'YLim',imgBigYLim)
    line([bLeftBig+mgScale bLeftBig+mgScale+scaleBar*1000/pixSize],...
        [bBottomBig+mgScale bBottomBig+mgScale],'LineWidth',2,'Color',[1,1,1]);
    text(bLeftBig+5, bBottomBig+mgScale+11,'1 \mum','Color',[1,1,1]);
    txtAx6= text(bRightBig-40, bBottomBig+12,[num2str(midEndFrameAdjusted*tInterval) ' sec'],'Color',[1,1,1]);
    set(txtAx6,'horizontalAlignment','right')
    set(txtAx6,'position',[bRightBig-mgScale,bBottomBig+12])
    %Montage
    cropImg = imgMap(bBottom:bTop,bLeft:bRight,chosenStartFrame:chosenEndFrame);
    cropTmap = tMap(bBottom:bTop,bLeft:bRight,chosenStartFrame:chosenEndFrame);
    % make 4D array to use montage
    cropImg = reshape(cropImg,size(cropImg,1),size(cropImg,2),1,size(cropImg,3));
    cropTmap = reshape(cropTmap,size(cropTmap,1),size(cropTmap,2),1,size(cropTmap,3));
    axes('Position',[marginX, 430/figHeight+marginY, 540/figWidth-marginX,50/figHeight]);
    % Fix the individual montage size and determine the number
    montWidth = 25; % in pixel
    maxMontageNum = floor(535/montWidth)*2;
    numChosenFrames = length(chosenStartFrame:chosenEndFrame);
    montInterval = ceil(numChosenFrames/maxMontageNum);
    hm1=montage(cropImg,'Size',[2, NaN],'DisplayRange',[], 'Indices', 1:montInterval:numChosenFrames);
    line([2 2+scaleBar*500/pixSize],...
        [3 3],'LineWidth',2,'Color',[1,1,1]);
    monImgH = floor(hm1.YData(2)/2);
    numCols = ceil(length(1:montInterval:numChosenFrames)/2);
    monImgW = floor(hm1.XData(2)/numCols);
    p=0;
    for ii= 1:montInterval:numChosenFrames
        p=p+1;
        iCol=mod(p-1,numCols);
        q=floor((p-1)/numCols);
        txtMont1= text(1+(iCol)*monImgW, (q+1)*monImgH-1,[num2str((ii-1)*tInterval) ' s'],'Color',[1,1,1]);
        set(txtMont1,'Fontsize',6)
        set(txtMont1,'horizontalAlignment','right')
        set(txtMont1,'position',[(iCol+1)*monImgW-0.5, (q+1)*monImgH-1.5])
    end
    
    ax7=subplot('Position',[marginX, 200/figHeight, 540/figWidth-marginX,50/figHeight]);
    tCropMax = max(curTrack.forceMag(chosenStartFrame:chosenEndFrame))*1.1;
    tCropMin = min(curTrack.forceMag(chosenStartFrame:chosenEndFrame));
    montage(cropTmap,'DisplayRange',[tCropMin, tCropMax],'Size',[2, NaN], 'Indices', 1:montInterval:numChosenFrames)
    p=0;
    for ii= 1:montInterval:numChosenFrames
        p=p+1;
        iCol=mod(p-1,numCols);
        q=floor((p-1)/numCols);
        txtMont2= text(1+(iCol)*monImgW, (q+1)*monImgH-1,[num2str((ii-1)*tInterval) ' s'],'Color','w');
        set(txtMont2,'Fontsize',6)
        set(txtMont2,'horizontalAlignment','right')
        set(txtMont2,'position',[(iCol+1)*monImgW-0.5, (q+1)*monImgH-1.5])
    end
    colormap(ax4,'jet')
    colormap(ax5,'jet')
    colormap(ax6,'jet')
    colormap(ax7,'jet')
    % intensity plot
%     subplot('Position',[0,0,0.5,0.33]), plot(chosenStartFrame:chosenEndFrame,curTrack.ampTotal(chosenStartFrame:chosenEndFrame))
%     subplot('Position',[0.5,0,0.5,0.33]), plot(chosenStartFrame:chosenEndFrame,curTrack.forceMag(chosenStartFrame:chosenEndFrame),'r')
    % ax 8: intensity time series
    ax8=axes('Position',[50/figWidth, 40/figHeight, 180/figWidth-marginX,140/figHeight]);
    plot((chosenStartFrame:chosenEndFrame)*tInterval,curTrack.ampTotal(chosenStartFrame:chosenEndFrame))
    xlabel('Time (s)'); ylabel('Fluorescence intensity (a.u.)')
    set(ax8,'FontSize',8)
    % force time series
    ax9=axes('Position',[300/figWidth, 40/figHeight, 180/figWidth-marginX,140/figHeight]);
    plot((chosenStartFrame:chosenEndFrame)*tInterval,curTrack.forceMag(chosenStartFrame:chosenEndFrame),'r')
    xlabel('Time (s)'); ylabel('Traction (Pa)')
    set(ax9,'FontSize',8)
    % force map colorbar
    ax10=axes('Position',[3*marginX+3*175/figWidth, 260/figHeight+marginY, 60/figWidth,145/figHeight]);
    axis tight
    caxis([tmin tmax]), axis off
%     if isempty(hc)
    hc = colorbar('West');
    set(hc,'FontSize',8)
    % force colorbar mini for montage
    ax11=axes('Position',[540/figWidth, 203/figHeight, 60/figWidth, 50/figHeight]);
    axis tight
    caxis([tCropMin, tCropMax]), axis off
%     if isempty(hc)
    hcMin = colorbar('West');
    set(hcMin,'FontSize',8)
    colormap(ax10,'jet')
    colormap(ax11,'jet')
    % saving
    iGroup = input('Which group does this adhesion belong to (1: turn-over 2: maturing 3: edge-transient)? ');
    if iGroup==1
        gPath = [outputPath filesep 'group1'];
        if ~exist(gPath,'dir')
            mkdir(gPath)
        end
    elseif iGroup==2
        gPath = [outputPath filesep 'group2'];
        if ~exist(gPath,'dir')
            mkdir(gPath)
        end
    elseif iGroup==3
        gPath = [outputPath filesep 'group3'];
        if ~exist(gPath,'dir')
            mkdir(gPath)
        end
    else
        disp('Unidentified group. Assining to group 4 ...')
        gPath = [outputPath filesep 'group4'];
        if ~exist(gPath,'dir')
            mkdir(gPath)
        end
    end
    print(h2,strcat(gPath,'/track',num2str(IDtoInspect),'.eps'),'-depsc2')
    savefig(h2,strcat(gPath,'/track',num2str(IDtoInspect),'.fig'))
    save([outputPath filesep 'selectedIDs.mat'], 'IDs')
    setappdata(hFig,'IDs',IDs);
    close(h2)
    axes(handles.axes1)
end
function txt=myupdateDC(~,event_obj)
    tracksNA = getappdata(hFig,'tracksNA');
    curPos = event_obj.Position;
    CurrentFrame = round((get(handles.SliderFrame,'Value')));
    idCurrent = arrayfun(@(x) logical(x.presence(CurrentFrame)),tracksNA);
    indexCur = find(idCurrent);
    idx = KDTreeClosestPoint([arrayfun(@(x) x.xCoord(CurrentFrame),tracksNA(idCurrent)),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNA(idCurrent))],curPos);
    selectedID = indexCur(idx);
    set(handles.Edit2,'String',num2str(selectedID));
    
    setappdata(hFig,'selPointID',selectedID); 
    try
        txt = {['ID: ', num2str(selectedID)],['Amp: ' num2str(tracksNA(selectedID).amp(CurrentFrame))],['Advance: ' num2str(tracksNA(selectedID).advanceDist)]};
    catch
        txt = {['ID: ', num2str(selectedID)],['Amp: ' num2str(tracksNA(selectedID).amp(CurrentFrame))]};
    end
end
function XListenerCallBack

    %// Retrieve handles structure. Used to let MATLAB recognize the
    %// edit box, slider and all UI components.
    handles = guidata(gcf);

    %// Here retrieve MyMatrix using getappdata.
    imgMap = getappdata(hFig,'MyMatrix');
    tracksNA = getappdata(hFig,'tracksNA');

    %// Get current frame
    CurrentFrame = round((get(handles.SliderFrame,'Value')));
    set(handles.Edit1,'String',num2str(CurrentFrame));
    prevXLim = handles.axes1.XLim;
    prevYLim = handles.axes1.YLim;

    %// Display appropriate frame.
    imshow(imgMap(:,:,CurrentFrame),[],'Parent',handles.axes1); 
    set(handles.axes1,'XLim',prevXLim,'YLim',prevYLim)
    hold on
    idCurrent = arrayfun(@(x) logical(x.presence(CurrentFrame)),tracksNA);
    plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNA(idCurrent)),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNA(idCurrent)),'ro')

    try
        idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{CurrentFrame}),tracksNA(idAdhLogic));
        idAdh = find(idAdhLogic);
        idAdhCur = idAdh(idAdhCur);
        arrayfun(@(x) plot(x.adhBoundary{CurrentFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idAdhCur))
    catch
        disp(' ')
    end
    
%     arrayfun(@(x) plot(x.adhBoundary{CurrentFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idCurrent))
    xmat = cell2mat(arrayfun(@(x) x.xCoord(1:CurrentFrame),tracksNA(idCurrent),'UniformOutput',false));
    ymat = cell2mat(arrayfun(@(x) x.yCoord(1:CurrentFrame),tracksNA(idCurrent),'UniformOutput',false));
    if size(xmat,2)==1
        plot(xmat',ymat','r.')
    else
        plot(xmat',ymat','r')
    end
    hold off

    guidata(hFig,handles);
end
%// Slider callback; executed when the slider is release or you press
%// the arrows.
function XSliderCallback(~,~)

    handles = guidata(gcf);

    %// Here retrieve MyMatrix using getappdata.
    imgMap = getappdata(hFig,'MyMatrix');
    tracksNA = getappdata(hFig,'tracksNA');

    CurrentFrame = round((get(handles.SliderFrame,'Value')));
    set(handles.Edit1,'String',num2str(CurrentFrame));
    prevXLim = handles.axes1.XLim;
    prevYLim = handles.axes1.YLim;

    imshow(imgMap(:,:,CurrentFrame),[],'Parent',handles.axes1); 
    set(handles.axes1,'XLim',prevXLim,'YLim',prevYLim)
    hold on
    idCurrent = arrayfun(@(x) logical(x.presence(CurrentFrame)),tracksNA);
    plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNA(idCurrent)),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNA(idCurrent)),'ro')
    
    idAdhLogic = arrayfun(@(x) ~isempty(x.adhBoundary),tracksNA);
    try
        idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{CurrentFrame}),tracksNA(idAdhLogic));
        idAdh = find(idAdhLogic);
        idAdhCur = idAdh(idAdhCur);
        arrayfun(@(x) plot(x.adhBoundary{CurrentFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idAdhCur))
    catch
        disp(' ')
    end
        
    xmat = cell2mat(arrayfun(@(x) x.xCoord(1:CurrentFrame),tracksNA(idCurrent),'UniformOutput',false));
    ymat = cell2mat(arrayfun(@(x) x.yCoord(1:CurrentFrame),tracksNA(idCurrent),'UniformOutput',false));
    if size(xmat,2)==1
        plot(xmat',ymat','r.')
    else
        plot(xmat',ymat','r')
    end
    hold off

    guidata(hFig,handles);
end
function windlowClose(~,~)
    IDs=getappdata(hFig,'IDs'); 
    if isempty(IDs)
        disp('No track was selected by data cursor. No ID is returend...')
    end
end
end
% for ii=startFrame:endFrame
%     p=p+1;
%     % actual frame
%     curFrame = imgMap(:,:,ii);
%     imshow(curFrame,[]), hold on
%     if ischar(idList) && strcmp(idList,'all')
%         plot(arrayfun(@(x) x.xCoord(ii),tracksNA),arrayfun(@(x) x.yCoord(ii),tracksNA),'ro')
%         xmat = cell2mat(arrayfun(@(x) x.xCoord(1:ii),tracksNA,'UniformOutput',false));
%         ymat = cell2mat(arrayfun(@(x) x.yCoord(1:ii),tracksNA,'UniformOutput',false));
%         plot(xmat',ymat','r')
%     else
%         plot(arrayfun(@(x) x.xCoord(ii),tracksNA(idList)),arrayfun(@(x) x.yCoord(ii),tracksNA(idList)),'ro')
%         xmat = cell2mat(arrayfun(@(x) x.xCoord(1:ii),tracksNA(idList),'UniformOutput',false));
%         ymat = cell2mat(arrayfun(@(x) x.yCoord(1:ii),tracksNA(idList),'UniformOutput',false));
%         plot(xmat',ymat','r')
%     end
%     drawnow
%     F(p) = getframe;
%     hold off
% end
% close(h)
% end
    %// Listener callback, executed when you drag the slider.

function drawImage(ii)
    % actual frame
    curFrame = imgMap(:,:,ii);
    imshow(curFrame,[]), hold on
    if ischar(idList) && strcmp(idList,'all')
        plot(arrayfun(@(x) x.xCoord(ii),tracksNA),arrayfun(@(x) x.yCoord(ii),tracksNA),'ro')
        xmat = cell2mat(arrayfun(@(x) x.xCoord(1:ii),tracksNA,'UniformOutput',false));
        ymat = cell2mat(arrayfun(@(x) x.yCoord(1:ii),tracksNA,'UniformOutput',false));
        plot(xmat',ymat','r')
    else
        plot(arrayfun(@(x) x.xCoord(ii),tracksNA(idList)),arrayfun(@(x) x.yCoord(ii),tracksNA(idList)),'ro')
        xmat = cell2mat(arrayfun(@(x) x.xCoord(1:ii),tracksNA(idList),'UniformOutput',false));
        ymat = cell2mat(arrayfun(@(x) x.yCoord(1:ii),tracksNA(idList),'UniformOutput',false));
        plot(xmat',ymat','r')
    end
end
        

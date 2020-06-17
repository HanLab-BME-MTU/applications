function [IDs, iGroups, iPrevGroups, tracksNA]=pickAdhesionTracksInteractive(tracksNA, imgMap, varargin)
% this function reads from Colocalization folder and shows the specified
% tracks on selected Channel interactively.
%% input reading
ip =inputParser;
ip.addRequired('tracksNA',@isstruct)
ip.addOptional('imgMap',[],@(x) isnumeric(x)); % This is the master channle index.
ip.addParamValue('numChan',2,@isscalar); % selcted track ids
ip.addParamValue('movieData',[],@(x) isa(x,'MovieData')); % selcted track ids
ip.addParamValue('trainedData',[],@istable); % trained data
ip.addParamValue('iChan',2,@isscalar); % This is the master channle index.
ip.addParamValue('iChanSlave',[],@(x) (isscalar(x) | isempty(x))); % This is the slave channle index.
ip.addParamValue('tMap',[],@(x) (isnumeric(x))); % This is the master channle index.
ip.addParamValue('imgMap2',[],@(x) (isnumeric(x))); % This is the master channle index.
ip.addParamValue('idSelected',[],@(x) (isstruct(x))); % This is the master channle index.

ip.parse(tracksNA, imgMap, varargin{:});
% pathForColocalization=ip.Results.pathForColocalization;
tracksNA=ip.Results.tracksNA;
numChan=ip.Results.numChan;
T=ip.Results.trainedData;
MD=ip.Results.movieData;
iChan=ip.Results.iChan;
iChanSlave=ip.Results.iChanSlave;
imgMap=ip.Results.imgMap;
tMap=ip.Results.tMap;
imgMap2=ip.Results.imgMap2;
idSelected=ip.Results.idSelected;
%% Load processed data
% movieData to find out pixel size
if isempty(MD)
    coloPath = fileparts(pathForColocalization);
    MDPath = fileparts(coloPath);
    MDfilePath = [MDPath filesep 'movieData.mat'];
    MD = load(MDfilePath,'MD');
    MD = MD.MD;
end
% iChan=2;
nFrames = MD.nFrames_;
tic
if isempty(tracksNA)
    tracksNA = load([pathForColocalization filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
    tracksNA = tracksNA.tracksNA;
end
if isempty(imgMap) %|| isempty(tMap)
    disp('Loading raw files ...')

    if numChan==1
        imgMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
        imgMap = imgMap.tMap;
    elseif numChan==2
        try
            imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
            imgMap = imgMap.paxImgStack;
            tMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
            tMap = tMap.tMap;
        catch
            % This case SDC was not used and first img frame was used.
            paxImage=MD.getChannel(iChan).loadImage(1); 
            [hp,w] = size(paxImage);
            if isempty(iChanSlave)
                TFMpackage=MD.getPackage(MD.getPackageIndex('TFMPackage'));
                forceProc =TFMpackage.processes_{4};
            %     forceField = forceProc.loadChannelOutput;
                tMapCell = load(forceProc.outFilePaths_{2});
                tMapCell = tMapCell.tMap;
                tMap = zeros(hp,w,nFrames);
                for iii=1:nFrames
                    tMap(:,:,iii) = tMapCell{iii};
                end
            else
                tMap = zeros(hp,w,nFrames);
                for iii=1:nFrames
                    slaveImage=MD.getChannel(iChanSlave).loadImage(iii); 
                    tMap(:,:,iii) = slaveImage;
                end
            end
            imgMap = zeros(hp,w,nFrames);
            for iii=1:nFrames
                paxImage=MD.getChannel(iChan).loadImage(iii); 
                imgMap(:,:,iii) = paxImage;
            end
        end
    end
end
outputPath = [MD.getPath filesep 'FAPackage' filesep 'trackAnalysis'];
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end
numFrames = size(imgMap,3);
startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA)));
endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA)));
pixSize = MD.pixelSize_; % nm/pixel
tInterval = MD.timeInterval_; % time interval in sec
scaleBar = 1; %micron

trainerInitially = false;
allDataClass = [];
% if ~isempty(T)
%     trainerInitially = true;
%     trainedClassifier = trainClassifierNA(T);
%     [~,allData] = extractFeatureNA(tracksNA);
%     allDataClass = predict(trainedClassifier,allData);
% end % This takes too much time SH 7/7/2016
% toc
%% show interactive movie interface
curNumFrames = endFrame-startFrame+1;
hFig = figure('Position',[100 100 720 750],'Units','normalized','DeleteFcn',@windlowClose);

handles.axes1 = axes('Units','normalized','Position',[0 0.05 1 0.95]);

%// Create slider and listener object for smooth visualization
handles.SliderFrame = uicontrol('Style','slider','Units','normalized','Position',[0 0 0.5 0.05],...
    'Min',startFrame,'Max',endFrame,'Value',startFrame,'SliderStep',[1/curNumFrames 2/curNumFrames],'Callback',@XSliderCallback);
handles.SliderxListener = handles.SliderFrame.addlistener('Value','PreSet',@(s,e) XListenerCallBack);
% handles.SliderxListener = addlistener(handles.SliderFrame,'Value','ContinuousValueChange',@(s,e) XListenerCallBack);

% handles.Text1 = uicontrol('Style','Text','Position',[180 420 60 30],'String','Current frame');
handles.Edit1 = uicontrol('Style','Edit','Units','normalized','Position',[0.5 0 0.05 0.05],'String',num2str(startFrame));
handles.Text2 = uicontrol('Style','Text','Units','normalized','Position',[0.55 0 0.15 0.05],'String',{'Adhesion ID ';'to inspect'});
handles.Edit2 = uicontrol('Style','Edit','Units','normalized','Position',[0.7 0 0.06 0.05],'String',num2str(startFrame));
handles.PushB2 = uicontrol('Style','pushbutton','Units','normalized','Position',[0.76 0 0.10 0.05],'String','Inspect','Callback',@pushInspectAdhesion);
handles.PushB3 = uicontrol('Style','pushbutton','Units','normalized','Position',[0.86 0 0.14 0.05],'String','RateConst','Callback',@pushRateConstant);

%// Use setappdata to store the image stack and in callbacks, use getappdata to retrieve it and use it. Check the docs for the calling syntax.
setappdata(hFig,'MyMatrix',imgMap); %// You could use %//setappdata(0,'MyMatrix',MyMatrix) to store in the base workspace. 
setappdata(hFig,'tracksNA',tracksNA); 
%// Display 1st frame
imshow(imgMap(:,:,startFrame),[]), hold on
if ~isempty(idSelected)
    htrackG = drawSelectedTracks(tracksNA,idSelected,1,gca);
elseif trainerInitially
    drawClassifiedTracks(allDataClass,tracksNA,1,gca,true);
else
%     idAdhLogic = arrayfun(@(x) ~isempty(x.adhBoundary),tracksNA);
%     idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{startFrame}),tracksNA(idAdhLogic));
%     idAdh = find(idAdhLogic);
%     idAdhCur = idAdh(idAdhCur);
%     arrayfun(@(x) plot(x.adhBoundary{startFrame}(:,1),x.adhBoundary{startFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idAdhCur))
    plot(arrayfun(@(x) x.xCoord(startFrame),tracksNA),arrayfun(@(x) x.yCoord(startFrame),tracksNA),'ro')
    xmat = cell2mat(arrayfun(@(x) x.xCoord(1:startFrame),tracksNA,'UniformOutput',false));
    ymat = cell2mat(arrayfun(@(x) x.yCoord(1:startFrame),tracksNA,'UniformOutput',false));
    if size(xmat,2)==1
        plot(xmat',ymat','r.')
    else
        plot(xmat',ymat','r')
    end
end
if ~isempty(idSelected)
    classDescription={'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
        'G4:retracting','G5:stable at the edge','G6:noise or very transient','G7:adhesions at stalling edge','G8:strong stable adhesion', 'G9:weak stable adhesion inside'};
    lgdHandle=legend([htrackG{1} htrackG{2} htrackG{3} htrackG{4} htrackG{5} htrackG{6} htrackG{7} htrackG{8} htrackG{9}],...
        classDescription(~cellfun(@isempty,htrackG)'),'TextColor','k','Location','best');
    lgdHandle.Color='k'; lgdHandle.TextColor='w';
end
hold off
% Supporting data cursor mode to identify an ID of NA track of interest.
dcm_obj = datacursormode(hFig);
set(dcm_obj,'UpdateFcn',@myupdateDC)
imgWidth = size(imgMap,2);
imgHeight = size(imgMap,1);
selectedID = [];
IDs=[];
iGroups=[];
iPrevGroups = [];
%// IMPORTANT. Update handles structure.
guidata(hFig,handles);
waitfor(hFig)
if isempty(IDs)
    disp('No track was selected by data cursor. No ID is returned...')
else
    disp(['Selected track is ' num2str(IDs) '.'])
end

function pushRateConstant(~,~)
    IDtoInspect=get(handles.Edit2,'String');
    IDtoInspect = str2double(IDtoInspect);
    curTrack = tracksNA(IDtoInspect);

    h2 = showSingleAdhesionTrackSummaryRateConstFitting(MD,curTrack,imgMap,tMap,imgMap2, IDtoInspect, outputPath);
    
end

function pushInspectAdhesion(~,~)
    IDtoInspect=get(handles.Edit2,'String');
    IDtoInspect = str2double(IDtoInspect);
    % show intensity profile and ask user to pick the starting frame and
    % ending frame
    % read first: 
    reAssign=false;
    newlyAssign=true;
    if ismember((IDtoInspect),IDs)
        disp(['The id, ' num2str(IDtoInspect) ' has been already selected for group ' num2str(iGroups(IDs==(IDtoInspect))) '. Do you want to reassign the group for this adhesion?((0)/1) ']);
        reAssign=input(['The id, ' num2str(IDtoInspect) ' has been already selected for group ' num2str(iGroups(IDs==(IDtoInspect))) '. Do you want to reassign the group for this adhesion?((0)/1) ']);
        if isempty(reAssign); reAssign=0; end
        whereInIDs = IDs==(IDtoInspect);
        newlyAssign = false;
    end
    if newlyAssign || reAssign
        if newlyAssign
            IDs=[IDs (IDtoInspect)];
        end
        curTrack = tracksNA(IDtoInspect);

        h2 = showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,imgMap2, IDtoInspect);

        % Display the features of the curTrack
        
        % Display the threshold of the currently assigned class
        
        
        % saving
        % This is a temporary fix for matlab 2018a, after this version this
        % disp function should be removed.
        disp({'Which group does this adhesion belong to :';  
            'group 1 (forming and disassembling as edge protrude),';
            'group 2 (maturing adhesions),';
            'group 3 (moving along with protruding edge)'; 'group 4 (adhesions at the retracting edges),';
            'group 5 (strong stable adhesion at the edge)'; 'group 6 (noisy and transient),';
            'group 7 (adhesions experiencing stalling edge)';
            'group 8 (strong stable adhesion inside edge) and )'; 
            'group 9 (weak adhesion inside edge)?'})
        iCurGroup = input(['Which group does this adhesion belong to :  \ngroup 1 (forming and disassembling as edge protrude),\ngroup 2 (maturing adhesions),', ...
            '\ngroup 3 (moving along with protruding edge), \ngroup 4 (adhesions at the retracting edges),' ...
            '\ngroup 5 (strong stable adhesion at the edge)','\ngroup 6 (noisy and transient) and ...','\ngroup 7 (adhesions experiencing stalling edge)',...
            '\ngroup 8 (strong stable adhesion inside edge) \ngroup 9 (weak adhesion inside edge)?']);
        if iCurGroup==1
            gPath = [outputPath filesep 'group1'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        elseif iCurGroup==2
            gPath = [outputPath filesep 'group2'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        elseif iCurGroup==3
            gPath = [outputPath filesep 'group3'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        elseif iCurGroup==4
            gPath = [outputPath filesep 'group4'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        elseif iCurGroup==5
            gPath = [outputPath filesep 'group5'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        elseif iCurGroup==6
            gPath = [outputPath filesep 'group6'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        elseif iCurGroup==7
            gPath = [outputPath filesep 'group7'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        elseif iCurGroup==8
            gPath = [outputPath filesep 'group8'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        elseif iCurGroup==9
            gPath = [outputPath filesep 'group9'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        else
            disp('Unidentified group. Assining to group 6 ...')
            iCurGroup=6;
            gPath = [outputPath filesep 'group6'];
            if ~exist(gPath,'dir')
                mkdir(gPath)
            end
        end
        if newlyAssign
            iGroups=[iGroups iCurGroup];
            if ~isempty(idSelected)
                iCurPrevGroup = 0;
                for kk=1:numel(fieldnames(idSelected))
                    if isfield(idSelected,'idGroup1Selected')
                        memberName = ['idGroup' num2str(kk) 'Selected'];
                        if ismember(IDtoInspect,idSelected.(memberName))
                            iCurPrevGroup = kk;
                            break
                        end
                    else
                        memberName = ['idGroup' num2str(kk)];
                        if ismember(IDtoInspect,find(idSelected.(memberName)))
                            iCurPrevGroup = kk;
                            break
                        end
                    end
                end
                iPrevGroups = [iPrevGroups iCurPrevGroup];
            end
%             idGroupSelected = sortIDTracks((IDtoInspect),iCurGroup,true);
%             [curT] = extractFeatureNA(tracksNA,idGroupSelected);
%             T = [T; curT];
        end
        if reAssign
            iGroups(whereInIDs)=iCurGroup;
            idGroupSelected = sortIDTracks((IDtoInspect),iCurGroup,true);
%             [curT] = extractFeatureNA(tracksNA,idGroupSelected);
%             T(whereInIDs,:) =curT;
        end
        disp(['Currently labeled groups: ' num2str(iGroups)])
        disp(['Currently labeled IDs: ' num2str(IDs)])
        disp(['Just labeled ID: ' num2str((IDtoInspect))])
        if trainerInitially
            disp('Training the classifier ...')
            tic
            trainedClassifier = trainClassifierNA(T);
        %     [~,allData] = extractFeatureNA(tracksNA);
            allDataClass = predict(trainedClassifier,allData);
            toc
        end    
        try
            tracksNA(IDtoInspect) = curTrack;
        catch
            tracksNA(IDtoInspect).startingFrameExtraExtra = [];
            tracksNA(IDtoInspect).endingFrameExtraExtra = [];
            tracksNA(IDtoInspect) = curTrack;
        end
%         if trainerInitially
%             print(h2,strcat(gPath,'/track',num2str(IDtoInspect),'.eps'),'-depsc2')
%             savefig(h2,strcat(gPath,'/track',num2str(IDtoInspect),'.fig'))
%         end
%         save([outputPath filesep 'selectedIDs.mat'], 'IDs', 'iGroups')
        setappdata(hFig,'IDs',IDs);
        setappdata(hFig,'iGroups',iGroups);
        close(h2)
        axes(handles.axes1)
    end
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
        splineParam=0.1;
        d = tracksNA(selectedID).ampTotal;
        tRange = tracksNA(selectedID).iFrame;
        d(d==0)=NaN;
        warning('off','SPLINES:CHCKXYWP:NaNs')
        try
            sd_spline= csaps(tRange,d,splineParam);
        catch
            d = tracksNA(selectedID).amp;
            d(tracksNA(selectedID).startingFrameExtraExtra:tracksNA(selectedID).endingFrameExtraExtra) = ...
                tracksNA(selectedID).ampTotal(tracksNA(selectedID).startingFrameExtraExtra:tracksNA(selectedID).endingFrameExtraExtra);
            sd_spline= csaps(tRange,d,splineParam);
        end
        sd=ppval(sd_spline,tRange);
        sd(isnan(d))=NaN;
        %         sd(isnan(d)) = NaN;
        % Find the maximum
        [~,curFrameMaxAmp]=nanmax(sd);
        timeToMaxInten = curFrameMaxAmp-tracksNA(selectedID).startingFrameExtra;
        
        
        txt = {['ID: ', num2str(selectedID)],...
            ['decayingIntensity: ' num2str(nanmax(tracksNA(selectedID).ampTotal)-tracksNA(selectedID).ampTotal(tracksNA(selectedID).endingFrameExtra))],...
            ['edgeAdvance: ' num2str(tracksNA(selectedID).edgeAdvanceDist(tracksNA(selectedID).endingFrameExtra))],...
            ['Advance: ' num2str(tracksNA(selectedID).advanceDist(tracksNA(selectedID).endingFrameExtra))],...
            ['lifeTime: ' num2str(tracksNA(selectedID).lifeTime)],...
            ['meanIntensity: ' num2str(nanmean(tracksNA(selectedID).amp))] , ...
            ['distToEdgeFirst: ' num2str(tracksNA(selectedID).distToEdge(tracksNA(selectedID).startingFrameExtra))], ...
            ['startingIntensity: ' num2str(tracksNA(selectedID).ampTotal(tracksNA(selectedID).startingFrameExtra))],...
            ['distToEdgeChange: ' num2str(tracksNA(selectedID).distToEdgeChange)]...
            ['distToEdgeLastNAs: ' num2str(tracksNA(selectedID).distToEdge(tracksNA(selectedID).endingFrameExtra))]...
            ['edgeAdvanceDistFirstChange: ' num2str(tracksNA(selectedID).advanceDistChange2min(min(tracksNA(selectedID).startingFrameExtra+30,tracksNA(selectedID).endingFrameExtra)))],...
            ['edgeAdvanceDistLastChange: ' num2str(tracksNA(selectedID).advanceDistChange2min(tracksNA(selectedID).endingFrameExtra))],...
            ['maxEdgeAdvanceDistChange: ' num2str(tracksNA(selectedID).maxEdgeAdvanceDistChange)],...
            ['maxIntensity: ' num2str(nanmax(tracksNA(selectedID).ampTotal))],...
            ['timeToMaxInten: ' num2str(timeToMaxInten)],...
            ['edgeVariation: ' num2str(min(nanstd(tracksNA(selectedID).closestBdPointNaive(:,1)),nanstd(tracksNA(selectedID).closestBdPointNaive(:,2))))]};
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
    if ~isempty(idSelected)
        htrackG = drawSelectedTracks(tracksNA,idSelected,CurrentFrame,gca);
    elseif trainerInitially
        drawClassifiedTracks(allDataClass(idCurrent,:),tracksNA(idCurrent),CurrentFrame,gca,true);
    else
%         arrayfun(@(x) plot(x.adhBoundary{CurrentFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idCurrent))
%         plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNA(idCurrent)),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNA(idCurrent)),'ro')
        xmat = cell2mat(arrayfun(@(x) x.xCoord(1:CurrentFrame),tracksNA(idCurrent),'UniformOutput',false));
        ymat = cell2mat(arrayfun(@(x) x.yCoord(1:CurrentFrame),tracksNA(idCurrent),'UniformOutput',false));
        if size(xmat,2)==1
            plot(xmat',ymat','r.')
        else
            plot(xmat',ymat','r')
        end
    end
%     try
%         idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{CurrentFrame}),tracksNA(idAdhLogic));
%         idAdh = find(idAdhLogic);
%         idAdhCur = idAdh(idAdhCur);
%         arrayfun(@(x) plot(x.adhBoundary{CurrentFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idAdhCur))
%     catch
%         disp(' ')
%     end
    if ~isempty(idSelected)
        lgdHandle=legend([htrackG{1} htrackG{2} htrackG{3} htrackG{4} htrackG{5} htrackG{6} htrackG{7} htrackG{8} htrackG{9}],...
            classDescription(~cellfun(@isempty,htrackG)'),'TextColor','k','Location','best');
        lgdHandle.Color='k'; lgdHandle.TextColor='w';
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
    zoom reset
    set(handles.axes1,'XLim',prevXLim,'YLim',prevYLim)
    hold on
    idCurrent = arrayfun(@(x) logical(x.presence(CurrentFrame)),tracksNA);
    if ~isempty(idSelected)
        htrackG = drawSelectedTracks(tracksNA,idSelected,CurrentFrame,gca);
    elseif trainerInitially
        drawClassifiedTracks(allDataClass(idCurrent,:),tracksNA(idCurrent),CurrentFrame,gca,true);
    else
%         arrayfun(@(x) plot(x.adhBoundary{CurrentFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idCurrent))
%         plot(arrayfun(@(x) x.xCoord(CurrentFrame),tracksNA(idCurrent)),arrayfun(@(x) x.yCoord(CurrentFrame),tracksNA(idCurrent)),'ro')
        xmat = cell2mat(arrayfun(@(x) x.xCoord(1:CurrentFrame),tracksNA(idCurrent),'UniformOutput',false));
        ymat = cell2mat(arrayfun(@(x) x.yCoord(1:CurrentFrame),tracksNA(idCurrent),'UniformOutput',false));
        if size(xmat,2)==1
            plot(xmat',ymat','r.')
        else
            plot(xmat',ymat','r')
        end
    end
    if ~isempty(idSelected)
        lgdHandle=legend([htrackG{1} htrackG{2} htrackG{3} htrackG{4} htrackG{5} htrackG{6} htrackG{7} htrackG{8} htrackG{9}],...
            classDescription(~cellfun(@isempty,htrackG)'),'TextColor','k','Location','best');
        lgdHandle.Color='k'; lgdHandle.TextColor='w';
    end
    %% segmented focal adhesions
%     idAdhLogic = arrayfun(@(x) ~isempty(x.adhBoundary),tracksNA);
%     try
%         idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{CurrentFrame}),tracksNA(idAdhLogic));
%         idAdh = find(idAdhLogic);
%         idAdhCur = idAdh(idAdhCur);
%         arrayfun(@(x) plot(x.adhBoundary{CurrentFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), 'Color',[255/255 153/255 51/255], 'LineWidth', 0.5),tracksNA(idAdhCur))
%     catch
%         disp(' ')
%     end
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
        

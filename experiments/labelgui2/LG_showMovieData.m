function figureHandle = LG_showMovieData(idlist, dataProperties, currentTime, loadedFrames, maxSpots, figureHandleOrPos, tagPopupMenuH, colorMap)
%LG_SHOWMOVIEDATA shows the movieData for labelgui2
%
% SYNOPSIS: figureHandle = LG_showMovieData(idlist, dataProperties, currentTime, loadedFrames, maxSpots, figureHandleOrPos)
%
% INPUT idlist: idlist
%		dataProperties: dataProperties
%		currentTime: currentTime
%		loadedFrames: list of movie frames currently loaded into memory
%       maxSpots: maximum number of spots loaded into memory
%		figureHandleOrPosition: if empty, or figurePosition, a new window will be opened.
%       tagPopupMenuH: handle to contextmenu of tags.
%       colorMap: colorMap used to plot tags
%
% OUTPUT figureHandle: handle to movieDataFigure
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 12-May-2006
%
% called by: LG_showMovieData_callback
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=======================
%% READ DATA FROM INPUT
%=======================

% we need
% - frameSize, frameSizeMu, pixelSize
% - loadedFrames (got that from input)
% - currentRealTime, interval
% - Status (either ok or deleted or tracked from #sourceList
% - max # of spots (got that from input)
% - for every spot
%   - name, flag, snr, amp, pos +/- uncertainty


% frameSize, frameSizeMu, pixelSize
pixelSize = [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z];
frameSize = dataProperties.movieSize(1:3);
frameSizeMu = frameSize .* pixelSize([1,1,2]);

% loadedFrames. Only keep first and last
loadedFrames = loadedFrames([1,end]);

% current real time
currentRealTime = mean(dataProperties.frameTime(currentTime,:));

% sampling
sampling = dataProperties.timeLapse;

% status
isDeleted = 0;
isTracked = 0;
if isempty(idlist)
    status = 'No idlist loaded';
else
if isempty(idlist(currentTime).linklist)
    status = 'Frame Deleted';
    isDeleted = 1;
elseif isfield(idlist(currentTime).info,'sourceList')
    status = sprintf(['Tracked from ', ...
        repmat('%i ',1,length(idlist(currentTime).info.sourceList))],...
        idlist(currentTime).info.sourceList);
    isTracked = 1;
elseif isfield(idlist(currentTime).info,'totalQ_Pix') && ...
        ~isempty(idlist(currentTime).info.totalQ_Pix)
    % this is for the few idlisttracks that were produced before the
    % sourceList was saved
    status = 'Tracked';
    isTracked = 1;
else
    status = '';
end
end

% loop through tags that belong to spots in current idlist to make
% movieInfo


if ~isempty(idlist) && ~isDeleted
    linklist = idlist(currentTime).linklist;
    nSpots = max(linklist(:,2));
    spotInfo(1:nSpots) = struct('name','','SNR','','amp','','flag','',...
        'x','','sx','','y','','sy','','z','','sz','',...
        'xmu','','sxmu','','ymu','','symu','','zmu','','szmu','',...
        'tagNumber','');

    for iTag = 1:size(linklist,1)
        % first, check if good spot. Then fill spotInfo
        if linklist(iTag,2) > 0 && linklist(iTag,3) ~= 4
            iSpot = linklist(iTag,2);
            spotInfo(iSpot).name = idlist(1).stats.labelcolor{iTag};
            spotInfo(iSpot).amp = linklist(iTag,8);
            spotInfo(iSpot).SNR = linklist(iTag,8)/sqrt(linklist(iTag,12));

            % x y z
            fn={'xmu','ymu','zmu';'sxmu','symu','szmu'};
            for dim=1:3
                pix = pixelSize((dim==3)+1);
            spotInfo(iSpot).(fn{1,dim}) = linklist(iTag,8+dim);
            if isTracked
                if isfield(idlist(currentTime).info,'totalQ_Pix')
                spotInfo(iSpot).(fn{2,dim}) = sqrt(...
                    idlist(currentTime).info.totalQ_Pix((iTag-1)*3+dim,(iTag-1)*3+dim) *...
                    idlist(currentTime).linklist(iTag,12)) * pix;
                else
                    spotInfo(iSpot).(fn{2,dim}) = sqrt(...
                    idlist(currentTime).info.detectQ_Pix((iTag-1)*3+dim,(iTag-1)*3+dim) *...
                    idlist(currentTime).linklist(iTag,12)) * pix;
                end
            else
                spotInfo(iSpot).(fn{2,dim}) = sqrt(...
                    idlist(currentTime).info.detectQ_Pix((iTag-1)*3+dim,(iTag-1)*3+dim) *...
                    idlist(currentTime).linklist(iTag,12)) * pix;
            end
            spotInfo(iSpot).(fn{1,dim}(1:end-2)) = ...
                spotInfo(iSpot).(fn{1,dim}) / pix;
            spotInfo(iSpot).(fn{2,dim}(1:end-2)) = ...
                spotInfo(iSpot).(fn{2,dim}) / pix;
            end
        

            % set flags
            flagList = {};
            if any(linklist(iTag,3) == [2,5]);
                flagList{end+1} = 'Trk';
            end
            if any(linklist(iTag,3) == [3,5]);
                flagList{end+1} = 'Fus';
            end
            if any(linklist(iTag,5) == [2]);
                flagList{end+1} = 'SO';
            end
            if isempty(flagList)
                flagList = {'--'};
            end
            spotInfo(iSpot).flag = flagList;

            % remember tagNumber
            spotInfo(iSpot).tagNumber = iTag;

        end
    end % loop tags

    % preserve order of tags, not spots.
    % therefore, sort spotInfo according to tagNumber
    tagNumbers = cat(1,spotInfo.tagNumber);
    [dummy,sortIdx] = sort(tagNumbers);
    spotInfo = spotInfo(sortIdx);
else
    % no spotInfo, b/c deleted tag
    nSpots = 0;
    spotInfo = [];
end

%=============================


%=============================
%% CREATE FIGURE
%=============================

if any(length(figureHandleOrPos) == [0,4]) || ~ishandle(figureHandleOrPos)
    % launch new figure.
    figureHandle = dialog('WindowStyle','normal',...
        'HandleVisibility','callback','Visible','off',...
        'Units','characters','Name','MovieData');
    position = get(figureHandle,'Position');
    % use input position if available
    if length(figureHandleOrPos) == 4
        position = figureHandleOrPos;
    end
    % set height, width
    position(3) = 50.2;
    height = 9 + 4* maxSpots;
    position(4) = height;
    set(figureHandle,'Position',position);

    % add textfields
    th.frame = uicontrol(figureHandle,'Style','text',...
        'Units','Characters', 'Position', [2,height-2.5,46.2,1.7],...
        'FontName','Helvetica','HorizontalAlignment','left',...
        'FontWeight','normal',...
        'String','');
    th.time = uicontrol(figureHandle,'Style','text',...
        'Units','Characters', 'Position', [2,height-4.5,46.2,1.7],...
        'FontName','Helvetica','HorizontalAlignment','left',...
        'FontWeight','normal','TooltipString',...
        't: current time. dt: average time between frames. Loaded: list of movie frames in memory',...
        'String','');
    th.status = uicontrol(figureHandle,'Style','text',...
        'Units','Characters', 'Position', [2,height-6.5,46.2,1.7],...
        'FontName','Helvetica','HorizontalAlignment','left',...
        'FontWeight','normal',...
        'String','');


    for iSpot = maxSpots:-1:1
        % make just two lines - don't make anything special for the name
        % right now
        th.spotName(iSpot) = uicontrol(figureHandle,'Style','text',...
            'Units','Characters', 'Position', [2,height-(9.5 + (iSpot-1)*4),12,1.5],...
            'FontName','Helvetica','HorizontalAlignment','center',...
            'FontWeight','bold',...
            'String','');
        th.spotAmp(iSpot) = uicontrol(figureHandle,'Style','text',...
            'Units','Characters', 'Position', [16,height-(9.5 + (iSpot-1)*4),32.2,1.5],...
            'FontName','Helvetica','HorizontalAlignment','left',...
            'FontWeight','normal',...
            'String','');
        th.spotPos(iSpot) = uicontrol(figureHandle,'Style','text',...
            'Units','Characters', 'Position', [2,height-(11 + (iSpot-1)*4),46.2,1.5],...
            'FontName','Helvetica','HorizontalAlignment','left',...
            'FontWeight','normal',...
            'String','');
    end

    set(figureHandle,'Visible','on')

    % remember handles
    setappdata(figureHandle,'textHandles',th);
    
    % remember maxSpots
    setappdata(figureHandle,'maxSpots',maxSpots);

    % add context menu to figure handles
    windowHandles = LG_createTagPopupMenu(guidata(figureHandle),figureHandle);
    guidata(figureHandle, windowHandles);
    tagPopupMenuH = windowHandles.tagPopupMenuH;

else
    % assign figureHandle
    figureHandle = figureHandleOrPos;
    % read handles to textfields
    th = getappdata(figureHandle,'textHandles');
end
figure(figureHandle)

%=======================



%=======================
%% WRITE INFORMATION
%=======================

% first line
set(th.frame,'String',sprintf(...
    'Frame %ix%ix%i   Pixelsize xy %2.3f z %2.3f',...
    frameSize, pixelSize),'TooltipString',sprintf(...
    'Frame in mu: %3.2f x %3.2f x %3.2f',frameSizeMu));

% second line
set(th.time,'String',sprintf(...
    't %4.2f  dt %1.3f   Loaded %i:%i',...
    currentRealTime,sampling, loadedFrames));

% third line
set(th.status,'String',status);

% spot info. First, clear all. Then write
if isfield(th,'spotName')
set([th.spotName;th.spotAmp;th.spotPos],'String','',...
        'uicontextmenu',[],...
        'UserData',[],'Visible','off');
end
   

for iSpot = 1:nSpots
    % first line
    set(th.spotName(iSpot),'String',spotInfo(iSpot).name,...
        'uicontextmenu',tagPopupMenuH,...
        'UserData',spotInfo(iSpot).tagNumber,'BackgroundColor',[0.5,0.5,0.5],...
        'ForegroundColor',colorMap(spotInfo(iSpot).tagNumber,:),'Visible','on');
    set(th.spotAmp(iSpot),'String',sprintf(...
        [repmat('%s ',length(spotInfo(iSpot).flag)),...
        ' SNR %3.2f  Amp %3.3f'],...
        spotInfo(iSpot).flag{:},...
        spotInfo(iSpot).SNR,spotInfo(iSpot).amp),...
        'TooltipString',...
        'Flags: -- none; Trk found by tracker; Fus fusion; SO single occurence',...
        'uicontextmenu',tagPopupMenuH,...
        'UserData',spotInfo(iSpot).tagNumber,'Visible','on');
    % second line
    set(th.spotPos(iSpot),'String',sprintf(...
        ['X %3.1f',char(177),'%1.1f  Y %3.1f',char(177),...
        '%1.1f  Z %3.1f',char(177),'%2.1f'],...
        spotInfo(iSpot).x,spotInfo(iSpot).sx,...
        spotInfo(iSpot).y,spotInfo(iSpot).sy,...
        spotInfo(iSpot).z,spotInfo(iSpot).sz),...
        'TooltipString',sprintf(...
        ['Positions in mu: X %3.3f',char(177),'%1.3f  Y %3.3f',char(177),...
        '%1.3f  Z %3.3f',char(177),'%2.3f'],...
        spotInfo(iSpot).xmu,spotInfo(iSpot).sxmu,...
        spotInfo(iSpot).ymu,spotInfo(iSpot).symu,...
        spotInfo(iSpot).zmu,spotInfo(iSpot).szmu),...
        'uicontextmenu',tagPopupMenuH,...
        'UserData',spotInfo(iSpot).tagNumber,'Visible','on');
end


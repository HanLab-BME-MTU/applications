function label_showMovieData(hObject,eventdata,handles,isRefresh)
%shows data about movie

if nargin<4 |isempty(isRefresh)
    isRefresh = 0;
end

%get handles
labelguiH = handles.labelgui;
labelguiHandles = handles;
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
dataFigH = GetUserData(imgFigureH,'dataFigH');

if isempty(dataFigH)&isRefresh
    dataFigH = myMessageBox(0,'initializing...','CurrentMovieData'); %no offset yet
end

if ~isRefresh
    %messagebox has to be launched or closed
    
    %if ischecked - uncheck
    isChecked = get(hObject, 'Checked');
    if strcmp(isChecked,'on')
        set(hObject, 'Checked', 'off');
        %close all windows
        dataFigH = findall(0,'Tag','CurrentMovieData');
        close(dataFigH);
        %delete all the saved movieDataFigs
        labelPanelList = findall(0,'Tag','LabelPanel');
        for i = 1:length(labelPanelList)
            SetUserData(labelPanelList(i),[],1,'CurrentMovieData');
        end
        return
    end
    
    %make sure there is a movie loaded
    if isempty(imgFigureH)
        disp('no movie loaded')
        return
        %no movie loaded, don't check
    end
    
    %else: isNotChecked - check, launch
    set(hObject, 'Checked', 'on');
    
    %launch
    dataFigH = myMessageBox(0,'initializing...','CurrentMovieData'); %no offset yet
    
    
    
else
    %make sure there is a movie loaded
    if isempty(imgFigureH)
        disp('no movie loaded')
        return
        %no movie loaded, keep checked though
    end
end

%save handle
SetUserData(imgFigureH,dataFigH,1,'CurrentMovieData');

%collect all necessary data

%labelguiData
curr_time = get(labelguiHandles.slider3,'Value');
sliceMode = get(labelguiHandles.radiobutton2, 'Value');
if sliceMode
    curr_slice = get(labelguiHandles.slider4, 'Value');
else
    curr_slice = [];
end

%labelPanelData
idlist = GetUserData(imgFigureH,'idlist');
dataProperties = GetUserData(imgFigureH,'dataProperties');
PIXELSIZE_XY=dataProperties.PIXELSIZE_XY;
PIXELSIZE_Z=dataProperties.PIXELSIZE_Z;
frameTime = dataProperties.frameTime;

%calc Message data
if sliceMode
    currentTime = frameTime(curr_time,curr_slice);
else
    currentTime = mean(frameTime(curr_time,:));
end
if isfield(idlist(curr_time).info,'trackerMessage')
    trackData = 1;
    if ~isempty(idlist(curr_time).info.trackerMessage.source)
        tracked = 1;
        trackerMessage = idlist(curr_time).info.trackerMessage;
        tMLength = length(trackerMessage);
    else
        tracked = 0;
        tMLength = 0;
    end
else
    trackData = 0;
    tracked = 0;
    tMLength = 0;
end

%compose Message
%get starting line
nameStart = 1;
frameStart = 2;
tagStart = frameStart + 3;
trackStart = tagStart + size(idlist(curr_time).linklist,1)*3+3;
nextStart = trackStart + trackData*3+tracked*(tMLength-3); 

%movie name
msg{1} = dataProperties.name;

%frame info
msg{frameStart} = ['current time [sec] = ',num2str(currentTime)];
msg{frameStart + 1} = ['pixelsize = ',num2str(PIXELSIZE_XY),' ',num2str(PIXELSIZE_Z)];
msg(frameStart + 2) = {''};

%tag info: write coords /n snr
if ~isempty(idlist(curr_time).linklist)
    msg{tagStart} = 'tag / x / y / z // amp / snr';
    for i = 1:size(idlist(curr_time).linklist,1)
        %coords
        coordString = sprintf('  x: %3.2f  y: %3.2f  z:%3.2f',...
                (idlist(curr_time).linklist(i,9:11)));
        msg{tagStart+3*i-1} = [idlist(1).stats.labelcolor{i},coordString];
        
        %snr & amp & Q
        ampString = sprintf('amp: %1.2d  ',idlist(curr_time).linklist(i,8));
        
        if size(idlist(curr_time).linklist,2)>11
            
            snrString = sprintf('snr: %2.2f',idlist(curr_time).linklist(i,8)/sqrt(idlist(curr_time).linklist(i,12)));
            msg{tagStart+3*i} = [ampString,snrString];
        else
            msg{tagStart+3*i} = [ampString,' snr: N/A'];
        end
    end
else
    msg{tagStart + 1} = 'frame deleted';
end

%track info
if trackData
    if ~tracked
        msg{trackStart} = 'not tracked';
        msg(trackStart+1:trackStart+2) = {''};
    else
        for i = 1:length(trackerMessage.source)
            msg{trackStart + 3*(i-1)} = ['Source: ', trackerMessage.source{i}];
            msg{trackStart + 1 + 3*(i-1)} = ['Number of Iterations: ', trackerMessage.iterations{i}];
            msg{trackStart + 2 + 3*(i-1)} = trackerMessage.message{i};
        end
    end
else
    msg(trackStart:trackStart+2) = {''};
end


%update messageBox
myMessageBox(dataFigH,msg);

%save messageBoxH
SetUserData(imgFigureH,dataFigH,1);
        

function label_updateIntFigCB(h,eventdata,handles)
%callback for PushButton updateIntFig_PB to update plot of tag intensities in intFigure

%get handles
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
labelGuiH=findall(0,'Tag','labelgui');

%check whether there is an idlist loaded
if isempty(imgFigureH)
    h = errordlg('load movie and idlist first!');
    uiwait(h);
    return
else
    idlist = GetUserData(imgFigureH,'idlist');
    if isempty(idlist)
        h = errordlg('load idlist first!');
        uiwait(h);
        return
    end
end

%now build an intList from idlist.linklist
%1) get data
nTimePoints = size(idlist,2);
nColors = size(idlist(1).stats.labelcolor,1);
goodTime = [1:nTimePoints]';

%2) init intList
intList = zeros(nTimePoints,nColors);

%3) read linklist into intList
%move backwards in time to delete entries from goodTime
for t = nTimePoints:-1:1
    if isempty(idlist(t).linklist)
        goodTime(t)=[];
    else
        sortedLL=sortrows(idlist(t).linklist,4);
        intList(t,:) = sortedLL(:,8)';
    end
end

%update figure

%make figure active
intFigH = GetUserData(imgFigureH,'intFigH');
figure(intFigH);

%get colormap and update colorOrder
cMap=GetUserData(labelGuiH,'cMap');
cMapFact=size(cMap,1)/idlist(1).stats.maxColor;
labelColor=cMap(bsum2bvec(idlist(1).stats.maxColor-1)*cMapFact,:);
set(gca,'ColorOrder',labelColor);

%delete old lines
if ~isempty(handles)
    if isfield(handles,'spotIntLineH')
        for i = 1:length(handles.spotIntLineH)
            if ishandle(handles.spotIntLineH(i))
                delete(handles.spotIntLineH(i))
            end
        end
    else
        %no lines to delete, nothing to redraw
        return
    end
else
    %no lineH saved, nothing to redraw
    return
end

%plot new lines
spotIntLineH = plot(goodTime,intList(goodTime,:));

%update the black total intensity dots
frameIntH = findall(intFigH,'Type','line','LineStyle','none');
%read their intData
xData = [];
yData = [];
for i = 1:length(frameIntH)
    xData = [xData;get(frameIntH(i),'XData')'];
    yData = [yData;get(frameIntH(i),'YData')'];
end
[xData,sortIdx] = sort(xData);
yData = yData(sortIdx);

%delete old
delete(frameIntH);

%predefine frameIntH; if there's nothing to plot, the respective entry will
%be deleted (plot returns [])
frameIntH = zeros(2,1);

%now: write + for good, o for bad
tmpH = plot(goodTime,sum(intList(goodTime,:),2),'+k');
if isempty(tmpH)
    frameIntH(1) = [];
else
    frameIntH(1) = tmpH;
end
%badTime: all that has been plotted before, but not anymore
[badTime,badIdx] = setdiff(xData,goodTime);

tmpH = plot(badTime,yData(badIdx),'ok');
if isempty(tmpH)
    frameIntH(end) = [];
else
    frameIntH(end) = tmpH;
end

%set buttonDownFcn for the new lines
set(spotIntLineH,'ButtonDownFcn','label_gotoFrame_BDFCN');
set(frameIntH,'ButtonDownFcn','label_gotoFrame_BDFCN');

% show current time if necessary
% Check toggle status
showTime = get(handles.showTimeH(1),'Value');
% remove old line if necessary
if handles.showTimeH(2)
    % we don't want to delete 0!
    delete(handles.showTimeH(2))
end
if showTime
    % plot new line from yMin to yMax. Therefore, read currentTime and
    % yLimits of axes
    currentTime = get(findall(labelGuiH,'Tag','slider3'),'Value');
    yLims = get(gca,'YLim');
    % plot red line
    handles.showTimeH(2) = plot([currentTime,currentTime],yLims,'-r');
    % arrange timeLine so that it appears below everything
    c=get(gca,'Children');
    set(gca,'Children',c([2:end,1]));
    
    % set buttonDownFcn for timeLine
    set(handles.showTimeH(2),'ButtonDownFcn','label_gotoFrame_BDFCN');
else
    % remember that there should be no showTimeH(2), and don't set BDFCN
    handles.showTimeH(2) = 0;
end




%save new line data
handles.spotIntLineH = spotIntLineH;
guidata(h,handles);
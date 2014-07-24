function label_showDisplacementFigure(hObject,eventdata,handles,isUpdate)
%if exist, loads and shows displacementFigure

%check input
if nargin < 4 | isempty(isUpdate)
    isUpdate = 0;
end


if ~isUpdate
    %is item checked
    isChecked = get(hObject,'Checked');

    %check whether menuItem is checked or not only if original call
    if strcmp(isChecked,'on') %uncheck, close all windows
        set(hObject,'Checked','off');
        disFigH = findall(0,'Tag','disFig');
        close(disFigH);
        %delete all the saved disFigH
        labelPanelList = findall(0,'Tag','LabelPanel');
        for i = 1:length(labelPanelList)
            SetUserData(labelPanelList(i),[],1,'disFigH');
        end
        return

    else %check and open figure
        set(hObject,'Checked','on');
    end
end



%get data & figureHandles
labelguiH = handles.labelgui;
imgFigureH = GetUserData(handles.labelgui,'currentWindow');

if isempty(imgFigureH)
    h = errordlg('No movie/idlist loaded!')
    uiwait(h);
    set(hObject,'Checked','off');
    return
end

dataProperties = GetUserData(imgFigureH,'dataProperties');
idlist = GetUserData(imgFigureH,'idlist');

%load disFigureHandle from labelPanels
disFigH = GetUserData(imgFigureH,'disFigH');
if ~isempty(disFigH) & ishandle(disFigH)
    figure(disFigH);
    cla;
else
    %there is no figure around. Open a completely new figure
    disFigH = figure('Tag','disFig','Name',dataProperties.name);
    isUpdate = 0;
end

nTags = length(idlist(1).stats.labelcolor);
if nTags == 1
    h = errordlg('can not display distances - only one tag');
    uiwait(h);
    %store handles
    delete(disFigH);
    disFigH = [];
    SetUserData(imgFigureH,disFigH,1,'disFigH');
    return
end

%calc anaDat. We don't care about lastResult
anaDat = adgui_calc_anaDat(idlist,dataProperties,'currentIdlist');

%-------plot data

%if only two tags, we plot the second respective to the first
%if more tags, we plot all vs. the centroid


%get nTags & switch
nTags = anaDat(1).info.nTags;
legendText = {};

switch nTags
    case 1
        %only one tag. Makes no point to draw, but let's do it - we can,
        %after all
        [yData,yLabel,legendText,dummy,yStats] =...
            adgui_calcPlotData_tagDisplacementCorrected([],anaDat,'y',legendText,dataProperties,0,1);
        %set figure title
        title(['displacement of tag ',anaDat(1).info.labelColor{1}]);
    case 2
        %2 tags. calculate tag centered displacement
        [yData,yLabel,legendText,dummy,yStats] =...
            adgui_calcPlotData_tagDisplacementCorrected([],anaDat,'y',legendText,dataProperties,2,[2,1]);

        %set figure title
        title(['displacement of tag ',anaDat(1).info.labelColor{1}, ' relative to tag ',anaDat(1).info.labelColor{2}]);

    otherwise
        %more tags. calculate centroid centered displacement
        for i = 1:nTags
            [yData(:,i),yLabel,legendText,dummy,yStats(i)] =...
                adgui_calcPlotData_tagDisplacementCorrected([],anaDat,'y',legendText,dataProperties,1,i);
        end

        %set figure title
        title(['displacement of all tags relative to centroid']);

end

%calculate time
%calc data
timePoints = cat(1,anaDat.timePoint);
%timePointsBF = (timePoints(1:end-1)+timePoints(2:end))/2;
% instead of using timepoints between frames, we use the first of the two
% timepoints
timePointsBF = timePoints(1:end-1);

%calculate number of timepoints between frames
numTimePointsBF = diff(timePoints);

xData = timePointsBF;
xLabel = 'Timepoints between frames';
xStats = [];

legendText = [{'TimepointsBF'};legendText];

% adjust distance by dividing by numTimePointsBF - if one, remains the
% same, if two, is halved etc.
yData = yData./repmat(numTimePointsBF,1,size(yData,2));


%use the right colorOrder. From MatlabHelp:
%Note that if the axes NextPlot property is set to
%replace (the default), high-level functions like plot reset the ColorOrder
%property before determining the colors to use. If you want MATLAB to use a
%ColorOrder that is different from the default, set NextPlot to
%replacechildren. You can also specify your own default ColorOrder.

%get colormap
cMap=GetUserData(labelguiH,'cMap');
cMapFact=size(cMap,1)/idlist(1).stats.maxColor;

colorOrder = cMap([1:nTags]*cMapFact,:);

axH = gca;
set(axH,'NextPlot','add','ColorOrder',colorOrder);

%plot data without errorbars (for now)

lineH = plot(xData,yData,'-d');

%add axes label
xlabel(xLabel);
ylabel(yLabel);

%store axesHandle, too
axesH = gca;

%-------end plot data

%-------add functionality to plot

%if not update: we add the refresh button (same as intFigureUpdateButton)
if ~isUpdate
    uh = uicontrol('Style','pushbutton',...
        'Tag','updateDisFig_PB',...
        'Position',[5,5,120,23],...
        'Callback','label_showDisplacementFigure(gcbo,[],guidata(openfig(''labelgui'',''reuse'')),1)',...
        'String','update displacements');
    % add showTimeButton
    uh2 = uicontrol('Style','togglebutton',...
        'Tag','showTimeDisFig_PB',...
        'Position',[130,5,60,23],...
        'Callback','label_showDisplacementFigure(gcbo,[],guidata(openfig(''labelgui'',''reuse'')),1)',...
        'String','show time');
else
    % we are updating. Try to find showTimeButton
    childH = get(disFigH,'Children');
    showTimeIdx = strmatch('showTime',get(childH,'Tag'));
    % find button status
    showTime = get(childH(showTimeIdx),'Value');

    % draw line only if button active. There is no need to delete lines,
    % because there is a cla command above
    if showTime
        % plot new line from yMin to yMax. Therefore, read currentTime and
        % yLimits of axes
        currentTime = get(findall(labelguiH,'Tag','slider3'),'Value');
        yLims = get(gca,'YLim');
        % plot red line
        lineH(end+1,1) = plot([currentTime,currentTime],yLims,'-r');
        % arrange timeLine so that it appears below everything
        c=get(gca,'Children');
        set(gca,'Children',c([2:end,1]));

    end
end
%set buttonDownFcn for axes and lines
set([lineH;axesH],'ButtonDownFcn','label_gotoFrame_BDFCN');

%-------end add functionality

%-------store handles & finish up

%make labelgui visible
if ishandle(imgFigureH)
    figure(imgFigureH);
end
figure(labelguiH);

%store handles
SetUserData(imgFigureH,disFigH,1,'disFigH');

%-------end store handles & finish up